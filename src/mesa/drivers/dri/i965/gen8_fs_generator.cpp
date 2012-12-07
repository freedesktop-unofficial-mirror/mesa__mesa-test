/*
 * Copyright © 2010, 2011, 2012 Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice (including the next
 * paragraph) shall be included in all copies or substantial portions of the
 * Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/** @file gen8_fs_generate.cpp
 *
 * Code generation for Gen8+ hardware.
 */

extern "C" {
#include "main/macros.h"
#include "brw_context.h"
} /* extern "C" */

#include "brw_fs.h"
#include "brw_cfg.h"
#include "glsl/ir_print_visitor.h"

gen8_fs_generator::gen8_fs_generator(struct brw_context *brw,
                                     struct brw_wm_compile *c,
                                     struct gl_shader_program *shader_prog,
                                     struct gl_fragment_program *fp,
                                     bool dual_source_output)
   : gen8_generator(brw, shader_prog, &fp->Base, c), c(c), fp(fp),
     dual_source_output(dual_source_output)
{
   shader =
      shader_prog ? shader_prog->_LinkedShaders[MESA_SHADER_FRAGMENT] : NULL;
}

gen8_fs_generator::~gen8_fs_generator()
{
}

void
gen8_fs_generator::mark_surface_used(unsigned surf_index)
{
   assert(surf_index < BRW_MAX_SURFACES);

   c->prog_data.base.binding_table.size_bytes =
      MAX2(c->prog_data.base.binding_table.size_bytes, (surf_index + 1) * 4);
}

void
gen8_fs_generator::generate_fb_write(fs_inst *ir)
{
   if (fp->UsesKill) {
      gen8_instruction *mov =
         MOV(retype(brw_vec1_grf(1, 7), BRW_REGISTER_TYPE_UW),
             brw_flag_reg(0, 1));
      mov->set_mask_control(BRW_MASK_DISABLE);
   }

   if (ir->header_present) {
      gen8_instruction *mov =
         MOV_RAW(brw_message_reg(ir->base_mrf), brw_vec8_grf(0, 0));
      mov->set_exec_size(BRW_EXECUTE_16);
   }

   gen8_instruction *inst = next_inst(BRW_OPCODE_SENDC);

   inst->set_dst(retype(vec8(brw_null_reg()), BRW_REGISTER_TYPE_UW));

   /* headerless version, just submit color payload */
   inst->set_src0(brw_vec8_grf(GEN7_MRF_HACK_START + ir->base_mrf, 0));

   /* Set up the "Message Specific Control" bits for the Data Port Message
    * Descriptor.  These are documented in the "Render Target Write" message's
    * "Message Descriptor" documentation (vol5c.2).
    */
   uint32_t msg_type;
   /* Set the Message Type */
   if (this->dual_source_output)
      msg_type = BRW_DATAPORT_RENDER_TARGET_WRITE_SIMD8_DUAL_SOURCE_SUBSPAN01;
   else if (dispatch_width == 16)
      msg_type = BRW_DATAPORT_RENDER_TARGET_WRITE_SIMD16_SINGLE_SOURCE;
   else
      msg_type = BRW_DATAPORT_RENDER_TARGET_WRITE_SIMD8_SINGLE_SOURCE_SUBSPAN01;

   uint32_t msg_control = msg_type;

   /* "Last Render Target Select" must be set on all writes to the last of
    * the render targets (if using MRT), or always for a single RT scenario.
    */
   if (ir->target == c->key.nr_color_regions - 1)
      msg_control |= (1 << 4); /* Last Render Target Select */

   uint32_t surf_index =
      c->prog_data.binding_table.render_target_start + ir->target;

   inst->set_dp_message(GEN6_SFID_DATAPORT_RENDER_CACHE,
                        surf_index,
                        GEN6_DATAPORT_WRITE_MESSAGE_RENDER_TARGET_WRITE,
                        msg_control,
                        ir->mlen,
                        0,
                        ir->header_present,
                        ir->eot);

   mark_surface_used(surf_index);
}

/* Computes the integer pixel x,y values from the origin.
 *
 * This is the basis of gl_FragCoord computation, but is also used
 * pre-gen6 for computing the deltas from v0 for computing
 * interpolation.
 */
void
gen8_fs_generator::generate_pixel_xy(struct brw_reg dst, bool is_x)
{
   struct brw_reg g1_uw = retype(brw_vec1_grf(1, 0), BRW_REGISTER_TYPE_UW);
   struct brw_reg src;
   struct brw_reg deltas;

   if (is_x) {
      src = stride(suboffset(g1_uw, 4), 2, 4, 0);
      deltas = brw_imm_v(0x10101010);
   } else {
      src = stride(suboffset(g1_uw, 5), 2, 4, 0);
      deltas = brw_imm_v(0x11001100);
   }

   if (dispatch_width == 16) {
      dst = vec16(dst);
   }

   /* We do this 8 or 16-wide, but since the destination is UW we
    * don't do compression in the 16-wide case.
    */
   gen8_instruction *inst = ADD(dst, src, deltas);
   inst->set_qtr_control(BRW_COMPRESSION_NONE);
}

void
gen8_fs_generator::generate_linterp(fs_inst *inst,
                                    struct brw_reg dst,
                                    struct brw_reg *src)
{
   struct brw_reg delta_x = src[0];
   struct brw_reg delta_y = src[1];
   struct brw_reg interp = src[2];

   (void) delta_y;
   assert(delta_y.nr == delta_x.nr + 1);
   PLN(dst, interp, delta_x);
}

void
gen8_fs_generator::generate_tex(fs_inst *ir,
                                struct brw_reg dst,
                                struct brw_reg src)
{
   int msg_type = -1;
   int rlen = 4;
   uint32_t simd_mode = BRW_SAMPLER_SIMD_MODE_SIMD8;

   if (dispatch_width == 16)
      simd_mode = BRW_SAMPLER_SIMD_MODE_SIMD16;

   switch (ir->opcode) {
   case SHADER_OPCODE_TEX:
      if (ir->shadow_compare) {
         msg_type = GEN5_SAMPLER_MESSAGE_SAMPLE_COMPARE;
      } else {
         msg_type = GEN5_SAMPLER_MESSAGE_SAMPLE;
      }
      break;
   case FS_OPCODE_TXB:
      if (ir->shadow_compare) {
         msg_type = GEN5_SAMPLER_MESSAGE_SAMPLE_BIAS_COMPARE;
      } else {
         msg_type = GEN5_SAMPLER_MESSAGE_SAMPLE_BIAS;
      }
      break;
   case SHADER_OPCODE_TXL:
      if (ir->shadow_compare) {
         msg_type = GEN5_SAMPLER_MESSAGE_SAMPLE_LOD_COMPARE;
      } else {
         msg_type = GEN5_SAMPLER_MESSAGE_SAMPLE_LOD;
      }
      break;
   case SHADER_OPCODE_TXS:
      msg_type = GEN5_SAMPLER_MESSAGE_SAMPLE_RESINFO;
      break;
   case SHADER_OPCODE_TXD:
      if (ir->shadow_compare) {
         msg_type = HSW_SAMPLER_MESSAGE_SAMPLE_DERIV_COMPARE;
      } else {
         msg_type = GEN5_SAMPLER_MESSAGE_SAMPLE_DERIVS;
      }
      break;
   case SHADER_OPCODE_TXF:
      msg_type = GEN5_SAMPLER_MESSAGE_SAMPLE_LD;
      break;
   case SHADER_OPCODE_TXF_MS:
      msg_type = GEN7_SAMPLER_MESSAGE_SAMPLE_LD2DMS;
      break;
   case SHADER_OPCODE_LOD:
      msg_type = GEN5_SAMPLER_MESSAGE_LOD;
      break;
   default:
      assert(!"not reached");
      break;
   }
   assert(msg_type != -1);

   if (simd_mode == BRW_SAMPLER_SIMD_MODE_SIMD16) {
      rlen = 8;
      dst = vec16(dst);
   }

   /* Load the message header if present.  If there's a texture offset,
    * we need to set it up explicitly and load the offset bitfield.
    * Otherwise, we can use an implied move from g0 to the first message reg.
    */
   if (ir->texture_offset) {
      /* Explicitly set up the message header by copying g0 to the MRF. */
      MOV_RAW(retype(brw_message_reg(ir->base_mrf), BRW_REGISTER_TYPE_UD),
              retype(brw_vec8_grf(0, 0), BRW_REGISTER_TYPE_UD));

      /* Then set the offset bits in DWord 2. */
      MOV_RAW(retype(brw_vec1_reg(BRW_MESSAGE_REGISTER_FILE,
                                  ir->base_mrf, 2), BRW_REGISTER_TYPE_UD),
              brw_imm_ud(ir->texture_offset));
   } else {
      assert(!ir->header_present);
      src = brw_message_reg(ir->base_mrf);
   }

   uint32_t surf_index =
      c->prog_data.base.binding_table.texture_start + ir->sampler;

   gen8_instruction *inst = next_inst(BRW_OPCODE_SEND);
   inst->set_dst(dst);
   inst->set_src0(src);
   inst->set_sampler_message(surf_index,
                             ir->sampler,
                             msg_type,
                             rlen,
                             ir->mlen,
                             ir->header_present,
                             simd_mode);

   mark_surface_used(surf_index);
}


/* For OPCODE_DDX and OPCODE_DDY, per channel of output we've got input
 * looking like:
 *
 * arg0: ss0.tl ss0.tr ss0.bl ss0.br ss1.tl ss1.tr ss1.bl ss1.br
 *
 * and we're trying to produce:
 *
 *           DDX                     DDY
 * dst: (ss0.tr - ss0.tl)     (ss0.tl - ss0.bl)
 *      (ss0.tr - ss0.tl)     (ss0.tr - ss0.br)
 *      (ss0.br - ss0.bl)     (ss0.tl - ss0.bl)
 *      (ss0.br - ss0.bl)     (ss0.tr - ss0.br)
 *      (ss1.tr - ss1.tl)     (ss1.tl - ss1.bl)
 *      (ss1.tr - ss1.tl)     (ss1.tr - ss1.br)
 *      (ss1.br - ss1.bl)     (ss1.tl - ss1.bl)
 *      (ss1.br - ss1.bl)     (ss1.tr - ss1.br)
 *
 * and add another set of two more subspans if in 16-pixel dispatch mode.
 *
 * For DDX, it ends up being easy: width = 2, horiz=0 gets us the same result
 * for each pair, and vertstride = 2 jumps us 2 elements after processing a
 * pair. But for DDY, it's harder, as we want to produce the pairs swizzled
 * between each other.  We could probably do it like ddx and swizzle the right
 * order later, but bail for now and just produce
 * ((ss0.tl - ss0.bl)x4 (ss1.tl - ss1.bl)x4)
 */
void
gen8_fs_generator::generate_ddx(fs_inst *inst,
                                struct brw_reg dst,
                                struct brw_reg src)
{
   struct brw_reg src0 = brw_reg(src.file, src.nr, 1,
                                 BRW_REGISTER_TYPE_F,
                                 BRW_VERTICAL_STRIDE_2,
                                 BRW_WIDTH_2,
                                 BRW_HORIZONTAL_STRIDE_0,
                                 BRW_SWIZZLE_XYZW, WRITEMASK_XYZW);
   struct brw_reg src1 = brw_reg(src.file, src.nr, 0,
                                 BRW_REGISTER_TYPE_F,
                                 BRW_VERTICAL_STRIDE_2,
                                 BRW_WIDTH_2,
                                 BRW_HORIZONTAL_STRIDE_0,
                                 BRW_SWIZZLE_XYZW, WRITEMASK_XYZW);
   ADD(dst, src0, negate(src1));
}

/* The negate_value boolean is used to negate the derivative computation for
 * FBOs, since they place the origin at the upper left instead of the lower
 * left.
 */
void
gen8_fs_generator::generate_ddy(fs_inst *inst,
                                struct brw_reg dst,
                                struct brw_reg src,
                                bool negate_value)
{
   struct brw_reg src0 = brw_reg(src.file, src.nr, 0,
                                 BRW_REGISTER_TYPE_F,
                                 BRW_VERTICAL_STRIDE_4,
                                 BRW_WIDTH_4,
                                 BRW_HORIZONTAL_STRIDE_0,
                                 BRW_SWIZZLE_XYZW, WRITEMASK_XYZW);
   struct brw_reg src1 = brw_reg(src.file, src.nr, 2,
                                 BRW_REGISTER_TYPE_F,
                                 BRW_VERTICAL_STRIDE_4,
                                 BRW_WIDTH_4,
                                 BRW_HORIZONTAL_STRIDE_0,
                                 BRW_SWIZZLE_XYZW, WRITEMASK_XYZW);
   if (negate_value)
      ADD(dst, src1, negate(src0));
   else
      ADD(dst, src0, negate(src1));
}

void
gen8_fs_generator::generate_scratch_write(fs_inst *inst, struct brw_reg dst)
{
   assert(inst->mlen != 0);
   assert(!"TODO: Implement generate_scratch_write.");
}

void
gen8_fs_generator::generate_scratch_read(fs_inst *inst, struct brw_reg dst)
{
   assert(inst->mlen != 0);
   assert(!"TODO: Implement generate_scratch_read.");
}

void
gen8_fs_generator::generate_scratch_read_gen7(fs_inst *inst, struct brw_reg dst)
{
   assert(inst->mlen != 0);
   assert(!"TODO: Implement generate_scratch_read_gen7.");
}

void
gen8_fs_generator::generate_uniform_pull_constant_load(fs_inst *inst,
                                                       struct brw_reg dst,
                                                       struct brw_reg index,
                                                       struct brw_reg offset)
{
   assert(inst->mlen == 0);

   assert(index.file == BRW_IMMEDIATE_VALUE &&
          index.type == BRW_REGISTER_TYPE_UD);
   uint32_t surf_index = index.dw1.ud;

   assert(offset.file == BRW_GENERAL_REGISTER_FILE);
   /* Reference only the dword we need lest we anger validate_reg() with
    * reg.width > reg.execszie.
    */
   offset = brw_vec1_grf(offset.nr, 0);

   gen8_instruction *send = next_inst(BRW_OPCODE_SEND);
   send->set_mask_control(BRW_MASK_DISABLE);

   /* We use the SIMD4x2 mode because we want to end up with 4 constants in
    * the destination loaded consecutively from the same offset (which appears
    * in the first component, and the rest are ignored).
    */
   dst.width = BRW_WIDTH_4;
   send->set_dst(dst);
   send->set_src0(offset);
   send->set_sampler_message(surf_index,
                             0, /* The LD message ignore the sampler unit. */
                             GEN5_SAMPLER_MESSAGE_SAMPLE_LD,
                             1, /* rlen */
                             1, /* mlen */
                             false, /* no header */
                             BRW_SAMPLER_SIMD_MODE_SIMD4X2);

   mark_surface_used(surf_index);
}

void
gen8_fs_generator::generate_varying_pull_constant_load(fs_inst *ir,
                                                       struct brw_reg dst,
                                                       struct brw_reg index,
                                                       struct brw_reg offset)
{
   /* Varying-offset pull constant loads are treated as a normal expression on
    * gen7, so the fact that it's a send message is hidden at the IR level.
    */
   assert(!ir->header_present);
   assert(!ir->mlen);

   assert(index.file == BRW_IMMEDIATE_VALUE &&
          index.type == BRW_REGISTER_TYPE_UD);
   uint32_t surf_index = index.dw1.ud;

   uint32_t simd_mode, rlen, mlen;
   if (dispatch_width == 16) {
      mlen = 2;
      rlen = 8;
      simd_mode = BRW_SAMPLER_SIMD_MODE_SIMD16;
   } else {
      mlen = 1;
      rlen = 4;
      simd_mode = BRW_SAMPLER_SIMD_MODE_SIMD8;
   }

   gen8_instruction *send = next_inst(BRW_OPCODE_SEND);
   send->set_dst(dst);
   send->set_src0(offset);
   send->set_sampler_message(surf_index,
                             0, /* The LD message ignore the sampler unit. */
                             GEN5_SAMPLER_MESSAGE_SAMPLE_LD,
                             rlen, /* rlen */
                             mlen, /* mlen */
                             false, /* no header */
                             simd_mode);

   mark_surface_used(surf_index);
}

/**
 * Cause the current pixel/sample mask (from R1.7 bits 15:0) to be transferred
 * into the flags register (f0.0).
 *
 * Used only on Gen6 and above.
 */
void
gen8_fs_generator::generate_mov_dispatch_to_flags(fs_inst *ir)
{
   struct brw_reg flags = brw_flag_reg(0, ir->flag_subreg);
   struct brw_reg dispatch_mask =
      retype(brw_vec1_grf(1, 7), BRW_REGISTER_TYPE_UW);

   gen8_instruction *mov = MOV(flags, dispatch_mask);
   mov->set_mask_control(BRW_MASK_DISABLE);
}

void
gen8_fs_generator::generate_discard_jump(fs_inst *ir)
{
   /* This HALT will be patched up at FB write time to point UIP at the end of
    * the program, and at brw_uip_jip() JIP will be set to the end of the
    * current block (or the program).
    */
   discard_halt_patches.push_tail(new(mem_ctx) ip_record(nr_inst));

   HALT();
}

void
gen8_fs_generator::patch_discard_jumps_to_fb_writes()
{
   if (discard_halt_patches.is_empty())
      return;

   /* There is a somewhat strange undocumented requirement of using
    * HALT, according to the simulator.  If some channel has HALTed to
    * a particular UIP, then by the end of the program, every channel
    * must have HALTed to that UIP.  Furthermore, the tracking is a
    * stack, so you can't do the final halt of a UIP after starting
    * halting to a new UIP.
    *
    * Symptoms of not emitting this instruction on actual hardware
    * included GPU hangs and sparkly rendering on the piglit discard
    * tests.
    */
   gen8_instruction *last_halt = HALT();
   last_halt->set_uip(16);
   last_halt->set_jip(16);

   int ip = nr_inst;

   foreach_list(node, &discard_halt_patches) {
      ip_record *patch_ip = (ip_record *) node;
      gen8_instruction *patch = &store[patch_ip->ip];
      assert(patch->opcode() == BRW_OPCODE_HALT);

      /* HALT takes an instruction distance from the pre-incremented IP. */
      patch->set_uip((ip - patch_ip->ip) * 16);
   }

   this->discard_halt_patches.make_empty();
}

/**
 * Sets the first dword of a vgrf for simd4x2 uniform pull constant
 * sampler LD messages.
 *
 * We don't want to bake it into the send message's code generation because
 * that means we don't get a chance to schedule the instruction.
 */
void
gen8_fs_generator::generate_set_simd4x2_offset(fs_inst *ir,
                                               struct brw_reg dst,
                                               struct brw_reg src,
                                               struct brw_reg value)
{
   assert(value.file == BRW_IMMEDIATE_VALUE);
   MOV_RAW(retype(brw_vec1_reg(dst.file, dst.nr, 0), value.type), value);
}

void
gen8_fs_generator::generate_code(exec_list *instructions)
{
   int last_native_inst_offset = next_inst_offset;
   const char *last_annotation_string = NULL;
   const void *last_annotation_ir = NULL;

   if (unlikely(INTEL_DEBUG & DEBUG_WM)) {
      if (shader) {
         printf("Native code for fragment shader %d (%d-wide dispatch):\n",
                shader_prog->Name, dispatch_width);
      } else {
         printf("Native code for fragment program %d (%d-wide dispatch):\n",
                prog->Id, dispatch_width);
      }
   }

   cfg_t *cfg = NULL;
   if (unlikely(INTEL_DEBUG & DEBUG_WM))
      cfg = new(mem_ctx) cfg_t(mem_ctx, instructions);

   foreach_list(node, instructions) {
      fs_inst *ir = (fs_inst *) node;
      struct brw_reg src[3], dst;

      if (unlikely(INTEL_DEBUG & DEBUG_WM)) {
         foreach_list(node, &cfg->block_list) {
            bblock_link *link = (bblock_link *)node;
            bblock_t *block = link->block;

            if (block->start == ir) {
               printf("   START B%d", block->block_num);
               foreach_list(predecessor_node, &block->parents) {
                  bblock_link *predecessor_link =
                     (bblock_link *)predecessor_node;
                  bblock_t *predecessor_block = predecessor_link->block;
                  printf(" <-B%d", predecessor_block->block_num);
               }
               printf("\n");
            }
         }

         if (last_annotation_ir != ir->ir) {
            last_annotation_ir = ir->ir;
            if (last_annotation_ir) {
               printf("   ");
               if (shader)
                  ((ir_instruction *) ir->ir)->print();
               else {
                  const prog_instruction *fpi;
                  fpi = (const prog_instruction *) ir->ir;
                  printf("%d: ", (int)(fpi - prog->Instructions));
                  _mesa_fprint_instruction_opt(stdout,
                                               fpi,
                                               0, PROG_PRINT_DEBUG, NULL);
               }
               printf("\n");
            }
         }
         if (last_annotation_string != ir->annotation) {
            last_annotation_string = ir->annotation;
            if (last_annotation_string)
               printf("   %s\n", last_annotation_string);
         }
      }

      for (unsigned int i = 0; i < 3; i++) {
         src[i] = brw_reg_from_fs_reg(&ir->src[i]);

         /* The accumulator result appears to get used for the
          * conditional modifier generation.  When negating a UD
          * value, there is a 33rd bit generated for the sign in the
          * accumulator value, so now you can't check, for example,
          * equality with a 32-bit value.  See piglit fs-op-neg-uvec4.
          */
         assert(!ir->conditional_mod ||
                ir->src[i].type != BRW_REGISTER_TYPE_UD ||
                !ir->src[i].negate);
      }
      dst = brw_reg_from_fs_reg(&ir->dst);

      default_state.conditional_mod = ir->conditional_mod;
      default_state.predicate = ir->predicate;
      default_state.predicate_inverse = ir->predicate_inverse;
      default_state.saturate = ir->saturate;
      default_state.flag_subreg_nr = ir->flag_subreg;

      if (dispatch_width == 16 && !ir->force_uncompressed)
         default_state.exec_size = BRW_EXECUTE_16;
      else
         default_state.exec_size = BRW_EXECUTE_8;

      /* fs_inst::force_sechalf is only used for original Gen4 code, so we
       * don't handle it.  Add qtr_control to default_state if that changes.
       */
      assert(!ir->force_sechalf);

      switch (ir->opcode) {
      case BRW_OPCODE_MOV:
         MOV(dst, src[0]);
         break;
      case BRW_OPCODE_ADD:
         ADD(dst, src[0], src[1]);
         break;
      case BRW_OPCODE_MUL:
         MUL(dst, src[0], src[1]);
         break;
      case BRW_OPCODE_MACH:
         MACH(dst, src[0], src[1]);
         break;

      case BRW_OPCODE_MAD:
         default_state.access_mode = BRW_ALIGN_16;
         if (dispatch_width == 16) {
            gen8_instruction *mad;
            mad = MAD(dst, src[0], src[1], src[2]);
            mad = MAD(sechalf(dst), sechalf(src[0]),
                      sechalf(src[1]), sechalf(src[2]));
            mad->set_qtr_control(GEN6_COMPRESSION_2Q);
         } else {
            MAD(dst, src[0], src[1], src[2]);
         }
         default_state.access_mode = BRW_ALIGN_1;
         break;

      case BRW_OPCODE_LRP:
         default_state.access_mode = BRW_ALIGN_16;
         if (dispatch_width == 16) {
            gen8_instruction *lrp;
            lrp = LRP(dst, src[0], src[1], src[2]);
            lrp = LRP(sechalf(dst), sechalf(src[0]),
                      sechalf(src[1]), sechalf(src[2]));
            lrp->set_qtr_control(GEN6_COMPRESSION_2Q);
         } else {
            LRP(dst, src[0], src[1], src[2]);
         }
         default_state.access_mode = BRW_ALIGN_1;
         break;


      case BRW_OPCODE_FRC:
         FRC(dst, src[0]);
         break;
      case BRW_OPCODE_RNDD:
         RNDD(dst, src[0]);
         break;
      case BRW_OPCODE_RNDE:
         RNDE(dst, src[0]);
         break;
      case BRW_OPCODE_RNDZ:
         RNDZ(dst, src[0]);
         break;

      case BRW_OPCODE_AND:
         AND(dst, src[0], src[1]);
         break;
      case BRW_OPCODE_OR:
         OR(dst, src[0], src[1]);
         break;
      case BRW_OPCODE_XOR:
         XOR(dst, src[0], src[1]);
         break;
      case BRW_OPCODE_NOT:
         NOT(dst, src[0]);
         break;
      case BRW_OPCODE_ASR:
         ASR(dst, src[0], src[1]);
         break;
      case BRW_OPCODE_SHR:
         SHR(dst, src[0], src[1]);
         break;
      case BRW_OPCODE_SHL:
         SHL(dst, src[0], src[1]);
         break;

      case BRW_OPCODE_F32TO16:
         F32TO16(dst, src[0]);
         break;
      case BRW_OPCODE_F16TO32:
         F16TO32(dst, src[0]);
         break;

      case BRW_OPCODE_CMP:
         CMP(dst, ir->conditional_mod, src[0], src[1]);
         break;
      case BRW_OPCODE_SEL:
         SEL(dst, src[0], src[1]);
         break;

      case BRW_OPCODE_BFREV:
         /* BFREV only supports UD type for src and dst. */
         BFREV(retype(dst, BRW_REGISTER_TYPE_UD),
               retype(src[0], BRW_REGISTER_TYPE_UD));
         break;

      case BRW_OPCODE_FBH:
         /* FBH only supports UD type for dst. */
         FBH(retype(dst, BRW_REGISTER_TYPE_UD), src[0]);
         break;

      case BRW_OPCODE_FBL:
         /* FBL only supports UD type for dst. */
         FBL(retype(dst, BRW_REGISTER_TYPE_UD), src[0]);
         break;

      case BRW_OPCODE_CBIT:
         /* CBIT only supports UD type for dst. */
         CBIT(retype(dst, BRW_REGISTER_TYPE_UD), src[0]);
         break;

      case BRW_OPCODE_ADDC: {
         gen8_instruction *addc = ADDC(dst, src[0], src[1]);
         addc->set_acc_wr_control(1);
         break;
      }

      case BRW_OPCODE_SUBB: {
         gen8_instruction *subb = SUBB(dst, src[0], src[1]);
         subb->set_acc_wr_control(1);
         break;
      }

      case BRW_OPCODE_BFE:
         default_state.access_mode = BRW_ALIGN_16;
         if (dispatch_width == 16) {
            gen8_instruction *bfe;
            bfe = BFE(dst, src[0], src[1], src[2]);
            bfe = BFE(sechalf(dst), sechalf(src[0]),
                      sechalf(src[1]), sechalf(src[2]));
            bfe->set_qtr_control(GEN6_COMPRESSION_2Q);
         } else {
            BFE(dst, src[0], src[1], src[2]);
         }
         default_state.access_mode = BRW_ALIGN_1;
         break;

      case BRW_OPCODE_BFI1:
         BFI1(dst, src[0], src[1]);
         break;

      case BRW_OPCODE_BFI2:
         default_state.access_mode = BRW_ALIGN_16;
         if (dispatch_width == 16) {
            gen8_instruction *bfi2;
            bfi2 = BFI2(dst, src[0], src[1], src[2]);
            bfi2 = BFI2(sechalf(dst), sechalf(src[0]),
                        sechalf(src[1]), sechalf(src[2]));
            bfi2->set_qtr_control(GEN6_COMPRESSION_2Q);
         } else {
            BFI2(dst, src[0], src[1], src[2]);
         }
         default_state.access_mode = BRW_ALIGN_1;
         break;

      case BRW_OPCODE_IF:
         IF(BRW_PREDICATE_NORMAL);
         break;

      case BRW_OPCODE_ELSE:
         ELSE();
         break;

      case BRW_OPCODE_ENDIF:
         ENDIF();
         break;

      case BRW_OPCODE_DO:
         DO();
         break;

      case BRW_OPCODE_BREAK:
         BREAK();
         break;

      case BRW_OPCODE_CONTINUE:
         CONTINUE();
         break;

      case BRW_OPCODE_WHILE:
         WHILE();
         break;

      case SHADER_OPCODE_RCP:
         MATH(BRW_MATH_FUNCTION_INV, dst, src[0]);
         break;

      case SHADER_OPCODE_RSQ:
         MATH(BRW_MATH_FUNCTION_RSQ, dst, src[0]);
         break;

      case SHADER_OPCODE_SQRT:
         MATH(BRW_MATH_FUNCTION_SQRT, dst, src[0]);
         break;

      case SHADER_OPCODE_EXP2:
         MATH(BRW_MATH_FUNCTION_EXP, dst, src[0]);
         break;

      case SHADER_OPCODE_LOG2:
         MATH(BRW_MATH_FUNCTION_LOG, dst, src[0]);
         break;

      case SHADER_OPCODE_SIN:
         MATH(BRW_MATH_FUNCTION_SIN, dst, src[0]);
         break;

      case SHADER_OPCODE_COS:
         MATH(BRW_MATH_FUNCTION_COS, dst, src[0]);
         break;

      case SHADER_OPCODE_INT_QUOTIENT:
         MATH(BRW_MATH_FUNCTION_INT_DIV_QUOTIENT, dst, src[0], src[1]);
         break;

      case SHADER_OPCODE_INT_REMAINDER:
         MATH(BRW_MATH_FUNCTION_INT_DIV_REMAINDER, dst, src[0], src[1]);
         break;

      case SHADER_OPCODE_POW:
         MATH(BRW_MATH_FUNCTION_POW, dst, src[0], src[1]);
         break;

      case FS_OPCODE_PIXEL_X:
         generate_pixel_xy(dst, true);
         break;
      case FS_OPCODE_PIXEL_Y:
         generate_pixel_xy(dst, false);
         break;
      case FS_OPCODE_CINTERP:
         MOV(dst, src[0]);
         break;
      case FS_OPCODE_LINTERP:
         generate_linterp(ir, dst, src);
         break;
      case SHADER_OPCODE_TEX:
      case FS_OPCODE_TXB:
      case SHADER_OPCODE_TXD:
      case SHADER_OPCODE_TXF:
      case SHADER_OPCODE_TXF_MS:
      case SHADER_OPCODE_TXL:
      case SHADER_OPCODE_TXS:
      case SHADER_OPCODE_LOD:
         generate_tex(ir, dst, src[0]);
         break;

      case SHADER_OPCODE_TG4:
      case SHADER_OPCODE_TG4_OFFSET:
         assert(!"XXX: Missing Gen8 scalar support for texture gather");
         break;

      case FS_OPCODE_DDX:
         generate_ddx(ir, dst, src[0]);
         break;
      case FS_OPCODE_DDY:
         /* Make sure fp->UsesDFdy flag got set (otherwise there's no
          * guarantee that c->key.render_to_fbo is set).
          */
         assert(fp->UsesDFdy);
         generate_ddy(ir, dst, src[0], c->key.render_to_fbo);
         break;

      case SHADER_OPCODE_GEN4_SCRATCH_WRITE:
         generate_scratch_write(ir, src[0]);
         break;

      case SHADER_OPCODE_GEN4_SCRATCH_READ:
         generate_scratch_read(ir, dst);
         break;

      case SHADER_OPCODE_GEN7_SCRATCH_READ:
         generate_scratch_read_gen7(ir, dst);
         break;

      case FS_OPCODE_UNIFORM_PULL_CONSTANT_LOAD_GEN7:
         generate_uniform_pull_constant_load(ir, dst, src[0], src[1]);
         break;

      case FS_OPCODE_VARYING_PULL_CONSTANT_LOAD_GEN7:
         generate_varying_pull_constant_load(ir, dst, src[0], src[1]);
         break;

      case FS_OPCODE_FB_WRITE:
         generate_fb_write(ir);
         break;

      case FS_OPCODE_MOV_DISPATCH_TO_FLAGS:
         generate_mov_dispatch_to_flags(ir);
         break;

      case FS_OPCODE_DISCARD_JUMP:
         generate_discard_jump(ir);
         break;

      case SHADER_OPCODE_SHADER_TIME_ADD:
         assert(!"XXX: Missing Gen8 scalar support for INTEL_DEBUG=shader_time");
         break;

      case SHADER_OPCODE_UNTYPED_ATOMIC:
         assert(!"XXX: Missing Gen8 scalar support for untyped atomics");
         break;

      case SHADER_OPCODE_UNTYPED_SURFACE_READ:
         assert(!"XXX: Missing Gen8 scalar support for untyped surface reads");
         break;

      case FS_OPCODE_SET_SIMD4X2_OFFSET:
         generate_set_simd4x2_offset(ir, dst, src[0], src[1]);
         break;

      case FS_OPCODE_SET_OMASK:
         assert(!"XXX: Missing Gen8 scalar support for SET_OMASK");
         break;

      case FS_OPCODE_SET_SAMPLE_ID:
         assert(!"XXX: Missing Gen8 scalar support for SET_SAMPLE_ID");
         break;

      case FS_OPCODE_PACK_HALF_2x16_SPLIT:
         assert(!"XXX: Missing Gen8 scalar support for PACK_HALF_2x16_SPLIT");
         break;

      case FS_OPCODE_UNPACK_HALF_2x16_SPLIT_X:
      case FS_OPCODE_UNPACK_HALF_2x16_SPLIT_Y:
         assert(!"XXX: Missing Gen8 scalar support for UNPACK_HALF_2x16_SPLIT");
         break;

      case FS_OPCODE_PLACEHOLDER_HALT:
         /* This is the place where the final HALT needs to be inserted if
          * we've emitted any discards.  If not, this will emit no code.
          */
         patch_discard_jumps_to_fb_writes();
         break;

      default:
         if (ir->opcode < int(ARRAY_SIZE(opcode_descs))) {
            _mesa_problem(ctx, "Unsupported opcode `%s' in FS",
                          opcode_descs[ir->opcode].name);
         } else {
            _mesa_problem(ctx, "Unsupported opcode %d in FS", ir->opcode);
         }
         abort();
      }

      if (unlikely(INTEL_DEBUG & DEBUG_WM)) {
         disassemble(stdout, last_native_inst_offset, next_inst_offset);

         foreach_list(node, &cfg->block_list) {
            bblock_link *link = (bblock_link *)node;
            bblock_t *block = link->block;

            if (block->end == ir) {
               printf("   END B%d", block->block_num);
               foreach_list(successor_node, &block->children) {
                  bblock_link *successor_link =
                     (bblock_link *)successor_node;
                  bblock_t *successor_block = successor_link->block;
                  printf(" ->B%d", successor_block->block_num);
               }
               printf("\n");
            }
         }
      }

      last_native_inst_offset = next_inst_offset;
   }

   if (unlikely(INTEL_DEBUG & DEBUG_WM)) {
      printf("\n");
   }

   patch_jump_targets();
}

const unsigned *
gen8_fs_generator::generate_assembly(exec_list *simd8_instructions,
                                     exec_list *simd16_instructions,
                                     unsigned *assembly_size)
{
   dispatch_width = 8;
   generate_code(simd8_instructions);

   if (simd16_instructions) {
      /* Align to a 64-byte boundary. */
      while ((nr_inst * sizeof(gen8_instruction)) % 64)
         NOP();

      /* Save off the start of this 16-wide program */
      c->prog_data.prog_offset_16 = nr_inst * sizeof(gen8_instruction);

      dispatch_width = 16;
      generate_code(simd16_instructions);
   }

   *assembly_size = next_inst_offset;
   return (const unsigned *) store;
}
