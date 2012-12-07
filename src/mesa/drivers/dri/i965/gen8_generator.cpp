/*
 * Copyright © 2012 Intel Corporation
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

/** @file gen8_generator.cpp
 *
 * Code generation for Gen8+ hardware, replacing the brw_eu_emit.c layer.
 */

extern "C" {
#include "main/compiler.h"
#include "main/macros.h"
#include "brw_context.h"
} /* extern "C" */

#include "glsl/ralloc.h"
#include "brw_eu.h"
#include "brw_reg.h"
#include "gen8_generator.h"

gen8_generator::gen8_generator(struct brw_context *brw,
                               struct gl_shader_program *shader_prog,
                               struct gl_program *prog,
                               void *mem_ctx)
   : shader_prog(shader_prog), prog(prog), brw(brw), mem_ctx(mem_ctx)
{
   ctx = &brw->ctx;

   memset(&default_state, 0, sizeof(default_state));
   default_state.mask_control = BRW_MASK_ENABLE;

   store_size = 1024;
   store = rzalloc_array(mem_ctx, gen8_instruction, store_size);
   nr_inst = 0;
   next_inst_offset = 0;

   /* Set up the control flow stacks. */
   if_stack_depth = 0;
   if_stack_array_size = 16;
   if_stack = rzalloc_array(mem_ctx, int, if_stack_array_size);

   loop_stack_depth = 0;
   loop_stack_array_size = 16;
   loop_stack = rzalloc_array(mem_ctx, int, loop_stack_array_size);
}

gen8_generator::~gen8_generator()
{
}

gen8_instruction *
gen8_generator::next_inst(unsigned opcode)
{
   gen8_instruction *inst;

   if (nr_inst + 1 > unsigned(store_size)) {
      store_size <<= 1;
      store = reralloc(mem_ctx, store, gen8_instruction, store_size);
      assert(store);
   }

   next_inst_offset += 16;
   inst = &store[nr_inst++];

   memset(inst, 0, sizeof(gen8_instruction));

   inst->set_opcode(opcode);
   inst->set_exec_size(default_state.exec_size);
   inst->set_access_mode(default_state.access_mode);
   inst->set_mask_control(default_state.mask_control);
   inst->set_cond_modifier(default_state.conditional_mod);
   inst->set_pred_control(default_state.predicate);
   inst->set_pred_inv(default_state.predicate_inverse);
   inst->set_saturate(default_state.saturate);
   inst->set_flag_subreg_nr(default_state.flag_subreg_nr);
   return inst;
}

#define ALU1(OP)                                           \
gen8_instruction *                                         \
gen8_generator::OP(struct brw_reg dst, struct brw_reg src) \
{                                                          \
   gen8_instruction *inst = next_inst(BRW_OPCODE_##OP);    \
   inst->set_dst(dst);                                     \
   inst->set_src0(src);                                    \
   return inst;                                            \
}

#define ALU2(OP)                                                             \
gen8_instruction *                                                           \
gen8_generator::OP(struct brw_reg dst, struct brw_reg s0, struct brw_reg s1) \
{                                                                            \
   gen8_instruction *inst = next_inst(BRW_OPCODE_##OP);                      \
   inst->set_dst(dst);                                                       \
   inst->set_src0(s0);                                                       \
   inst->set_src1(s1);                                                       \
   return inst;                                                              \
}

#define ALU3(OP)                                          \
gen8_instruction *                                        \
gen8_generator::OP(struct brw_reg dst, struct brw_reg s0, \
                   struct brw_reg s1,  struct brw_reg s2) \
{                                                         \
   return alu3(BRW_OPCODE_##OP, dst, s0, s1, s2);         \
}

#define ALU3F(OP) \
gen8_instruction *                                        \
gen8_generator::OP(struct brw_reg dst, struct brw_reg s0, \
                   struct brw_reg s1,  struct brw_reg s2) \
{                                                         \
   assert(dst.type == BRW_REGISTER_TYPE_F);               \
   assert(s0.type == BRW_REGISTER_TYPE_F);                \
   assert(s1.type == BRW_REGISTER_TYPE_F);                \
   assert(s2.type == BRW_REGISTER_TYPE_F);                \
   return alu3(BRW_OPCODE_##OP, dst, s0, s1, s2);         \
}

ALU2(ADD)
ALU2(AND)
ALU2(ASR)
ALU3(BFE)
ALU2(BFI1)
ALU3(BFI2)
ALU1(F32TO16)
ALU1(F16TO32)
ALU1(BFREV)
ALU1(CBIT)
ALU2(ADDC)
ALU2(SUBB)
ALU2(DP2)
ALU2(DP3)
ALU2(DP4)
ALU2(DPH)
ALU1(FBH)
ALU1(FBL)
ALU1(FRC)
ALU2(LINE)
ALU3F(LRP)
ALU3F(MAD)
ALU2(MUL)
ALU1(MOV)
ALU1(NOT)
ALU2(OR)
ALU2(PLN)
ALU1(RNDD)
ALU1(RNDE)
ALU1(RNDZ)
ALU2(SEL)
ALU2(SHL)
ALU2(SHR)
ALU2(XOR)

gen8_instruction *
gen8_generator::CMP(struct brw_reg dst, unsigned conditional,
                    struct brw_reg src0, struct brw_reg src1)
{
   gen8_instruction *inst = next_inst(BRW_OPCODE_CMP);
   inst->set_cond_modifier(conditional);
   inst->set_dst(dst);
   inst->set_src0(src0);
   inst->set_src1(src1);
   return inst;
}

gen8_instruction *
gen8_generator::MAC(struct brw_reg d, struct brw_reg s0, struct brw_reg s1)
{
   gen8_instruction *inst = next_inst(BRW_OPCODE_MAC);
   inst->set_dst(d);
   inst->set_src0(s0);
   inst->set_src1(s1);
   inst->set_acc_wr_control(true);
   return inst;
}

gen8_instruction *
gen8_generator::MACH(struct brw_reg d, struct brw_reg s0, struct brw_reg s1)
{
   gen8_instruction *inst = next_inst(BRW_OPCODE_MACH);
   inst->set_dst(d);
   inst->set_src0(s0);
   inst->set_src1(s1);
   inst->set_acc_wr_control(true);
   return inst;
}

static int
get_3src_subreg_nr(struct brw_reg reg)
{
   if (reg.vstride == BRW_VERTICAL_STRIDE_0) {
      assert(brw_is_single_value_swizzle(reg.dw1.bits.swizzle));
      return reg.subnr / 4 + BRW_GET_SWZ(reg.dw1.bits.swizzle, 0);
   } else {
      return reg.subnr / 4;
   }
}

gen8_instruction *
gen8_generator::alu3(unsigned opcode,
                     struct brw_reg dst,
                     struct brw_reg src0,
                     struct brw_reg src1,
                     struct brw_reg src2)
{
   /* MRFs haven't existed since Gen7, so we better not be using them. */
   if (dst.file == BRW_MESSAGE_REGISTER_FILE) {
      dst.file = BRW_GENERAL_REGISTER_FILE;
      dst.nr += GEN7_MRF_HACK_START;
   }

   gen8_instruction *inst = next_inst(opcode);
   assert(inst->access_mode() == BRW_ALIGN_16);

   assert(dst.file == BRW_GENERAL_REGISTER_FILE);
   assert(dst.nr < 128);
   assert(dst.address_mode == BRW_ADDRESS_DIRECT);
   assert(dst.type == BRW_REGISTER_TYPE_F ||
          dst.type == BRW_REGISTER_TYPE_D ||
          dst.type == BRW_REGISTER_TYPE_UD);
   inst->set_dst_3src_reg_nr(dst.nr);
   inst->set_dst_3src_subreg_nr(dst.subnr / 16);
   inst->set_dst_3src_writemask(dst.dw1.bits.writemask);
   inst->set_exec_size(BRW_EXECUTE_8);

   assert(src0.file == BRW_GENERAL_REGISTER_FILE);
   assert(src0.address_mode == BRW_ADDRESS_DIRECT);
   assert(src0.nr < 128);
   inst->set_src0_3src_swizzle(src0.dw1.bits.swizzle);
   inst->set_src0_3src_subreg_nr(get_3src_subreg_nr(src0));
   inst->set_src0_3src_rep_ctrl(src0.vstride == BRW_VERTICAL_STRIDE_0);
   inst->set_src0_3src_reg_nr(src0.nr);
   inst->set_src0_3src_abs(src0.abs);
   inst->set_src0_3src_negate(src0.negate);

   assert(src1.file == BRW_GENERAL_REGISTER_FILE);
   assert(src1.address_mode == BRW_ADDRESS_DIRECT);
   assert(src1.nr < 128);
   inst->set_src1_3src_swizzle(src1.dw1.bits.swizzle);
   inst->set_src1_3src_subreg_lo(get_3src_subreg_nr(src1) & 3);
   inst->set_src1_3src_subreg_hi(get_3src_subreg_nr(src1) >> 2);
   inst->set_src1_3src_rep_ctrl(src1.vstride == BRW_VERTICAL_STRIDE_0);
   inst->set_src1_3src_reg_nr(src1.nr);
   inst->set_src1_3src_abs(src1.abs);
   inst->set_src1_3src_negate(src1.negate);

   assert(src2.file == BRW_GENERAL_REGISTER_FILE);
   assert(src2.address_mode == BRW_ADDRESS_DIRECT);
   assert(src2.nr < 128);
   inst->set_src2_3src_swizzle(src2.dw1.bits.swizzle);
   inst->set_src2_3src_subreg_nr(get_3src_subreg_nr(src2));
   inst->set_src2_3src_rep_ctrl(src2.vstride == BRW_VERTICAL_STRIDE_0);
   inst->set_src2_3src_reg_nr(src2.nr);
   inst->set_src2_3src_abs(src2.abs);
   inst->set_src2_3src_negate(src2.negate);

   /* Set both the source and destination types based on dst.type, ignoring
    * the source register types.  The MAD and LRP emitters both ensure that
    * all register types are float.  The BFE and BFI2 emitters, however, may
    * send us mixed D and UD source types and want us to ignore that.
    */
   switch (dst.type) {
   case BRW_REGISTER_TYPE_F:
      inst->set_src_3src_type(BRW_3SRC_TYPE_F);
      inst->set_dst_3src_type(BRW_3SRC_TYPE_F);
      break;
   case BRW_REGISTER_TYPE_D:
      inst->set_src_3src_type(BRW_3SRC_TYPE_D);
      inst->set_dst_3src_type(BRW_3SRC_TYPE_D);
      break;
   case BRW_REGISTER_TYPE_UD:
      inst->set_src_3src_type(BRW_3SRC_TYPE_UD);
      inst->set_dst_3src_type(BRW_3SRC_TYPE_UD);
      break;
   }

   return inst;
}

gen8_instruction *
gen8_generator::math(unsigned math_function,
                     struct brw_reg dst,
                     struct brw_reg src0)
{
   gen8_instruction *inst = next_inst(BRW_OPCODE_MATH);

   assert(dst.file == BRW_GENERAL_REGISTER_FILE);
   assert(src0.file == BRW_GENERAL_REGISTER_FILE);
   assert(dst.hstride == BRW_HORIZONTAL_STRIDE_1);

   inst->set_math_function(math_function);
   inst->set_dst(dst);
   inst->set_src0(src0);
   return inst;
}

gen8_instruction *
gen8_generator::MATH(unsigned math_function,
                     struct brw_reg dst,
                     struct brw_reg src0)
{
   assert(src0.type == BRW_REGISTER_TYPE_F);
   gen8_instruction *inst = math(math_function, dst, src0);
   return inst;
}

gen8_instruction *
gen8_generator::MATH(unsigned math_function,
                     struct brw_reg dst,
                     struct brw_reg src0,
                     struct brw_reg src1)
{
   bool int_math =
      math_function == BRW_MATH_FUNCTION_INT_DIV_QUOTIENT ||
      math_function == BRW_MATH_FUNCTION_INT_DIV_REMAINDER ||
      math_function == BRW_MATH_FUNCTION_INT_DIV_QUOTIENT_AND_REMAINDER;

   if (int_math) {
      assert(src0.type != BRW_REGISTER_TYPE_F);
      assert(src1.type != BRW_REGISTER_TYPE_F);
   } else {
      assert(math_function == BRW_MATH_FUNCTION_POW);
      assert(src0.type == BRW_REGISTER_TYPE_F);
   }

   gen8_instruction *inst = math(math_function, dst, src0);
   inst->set_src1(src1);
   return inst;
}

gen8_instruction *
gen8_generator::MOV_RAW(struct brw_reg dst, struct brw_reg src0)
{
   gen8_instruction *inst = next_inst(BRW_OPCODE_MOV);
   inst->set_dst(retype(dst, BRW_REGISTER_TYPE_UD));
   inst->set_src0(retype(src0, BRW_REGISTER_TYPE_UD));
   inst->set_mask_control(BRW_MASK_DISABLE);

   return inst;
}


gen8_instruction *
gen8_generator::NOP()
{
   return next_inst(BRW_OPCODE_NOP);
}

void
gen8_generator::push_if_stack(gen8_instruction *inst)
{
   if_stack[if_stack_depth] = inst - store;

   ++if_stack_depth;
   if (if_stack_array_size <= if_stack_depth) {
      if_stack_array_size *= 2;
      if_stack = reralloc(mem_ctx, if_stack, int, if_stack_array_size);
   }
}

gen8_instruction *
gen8_generator::pop_if_stack()
{
   --if_stack_depth;
   return &store[if_stack[if_stack_depth]];
}

/**
 * Patch the IF and ELSE instructions to set the jump offsets (JIP and UIP.)
 */
void
gen8_generator::patch_IF_ELSE(gen8_instruction *if_inst,
                              gen8_instruction *else_inst,
                              gen8_instruction *endif_inst)
{
   assert(if_inst != NULL && if_inst->opcode() == BRW_OPCODE_IF);
   assert(else_inst == NULL || else_inst->opcode() == BRW_OPCODE_ELSE);
   assert(endif_inst != NULL && endif_inst->opcode() == BRW_OPCODE_ENDIF);

   endif_inst->set_exec_size(if_inst->exec_size());

   if (else_inst == NULL) {
      /* Patch IF -> ENDIF */
      if_inst->set_jip(16 * (endif_inst - if_inst));
      if_inst->set_uip(16 * (endif_inst - if_inst));
   } else {
      else_inst->set_exec_size(if_inst->exec_size());

      /* Patch IF -> ELSE and ELSE -> ENDIF:
       *
       * The IF's JIP should point at the instruction after the ELSE.
       * The IF's UIP should point to the ENDIF.
       *
       * Both are expressed in bytes, hence the multiply by 16...128-bits.
       */
      if_inst->set_jip(16 * (else_inst - if_inst + 1));
      if_inst->set_uip(16 * (endif_inst - if_inst));

      /* Patch ELSE -> ENDIF:
       *
       * Since we don't set branch_ctrl, both JIP and UIP point to ENDIF.
       */
      else_inst->set_jip(16 * (endif_inst - else_inst));
      else_inst->set_uip(16 * (endif_inst - else_inst));
   }
   endif_inst->set_jip(16);
}

gen8_instruction *
gen8_generator::IF(unsigned predicate)
{
   gen8_instruction *inst = next_inst(BRW_OPCODE_IF);
   inst->set_dst(vec1(retype(brw_null_reg(), BRW_REGISTER_TYPE_D)));
   inst->set_exec_size(default_state.exec_size);
   inst->set_pred_control(predicate);
   inst->set_mask_control(BRW_MASK_ENABLE);
   push_if_stack(inst);

   return inst;
}

gen8_instruction *
gen8_generator::ELSE()
{
   gen8_instruction *inst = next_inst(BRW_OPCODE_ELSE);
   inst->set_dst(retype(brw_null_reg(), BRW_REGISTER_TYPE_D));
   inst->set_mask_control(BRW_MASK_ENABLE);
   push_if_stack(inst);
   return inst;
}

gen8_instruction *
gen8_generator::ENDIF()
{
   gen8_instruction *if_inst = NULL;
   gen8_instruction *else_inst = NULL;

   gen8_instruction *tmp = pop_if_stack();
   if (tmp->opcode() == BRW_OPCODE_ELSE) {
      else_inst = tmp;
      tmp = pop_if_stack();
   }
   assert(tmp->opcode() == BRW_OPCODE_IF);
   if_inst = tmp;

   gen8_instruction *endif_inst = next_inst(BRW_OPCODE_ENDIF);
   endif_inst->set_mask_control(BRW_MASK_ENABLE);
   patch_IF_ELSE(if_inst, else_inst, endif_inst);

   return endif_inst;
}

unsigned
gen8_generator::next_ip(unsigned ip) const
{
   return ip + 16;
}

unsigned
gen8_generator::find_next_block_end(unsigned start) const
{
   for (unsigned ip = next_ip(start); ip < next_inst_offset; ip = next_ip(ip)) {
      gen8_instruction *inst = &store[ip / 16];

      switch (inst->opcode()) {
      case BRW_OPCODE_ENDIF:
      case BRW_OPCODE_ELSE:
      case BRW_OPCODE_WHILE:
      case BRW_OPCODE_HALT:
         return ip;
      }
   }

   return 0;
}

/* There is no DO instruction on Gen6+, so to find the end of the loop
 * we have to see if the loop is jumping back before our start
 * instruction.
 */
unsigned
gen8_generator::find_loop_end(unsigned start) const
{
   /* Always start after the instruction (such as a WHILE) we're trying to fix
    * up.
    */
   for (unsigned ip = next_ip(start); ip < next_inst_offset; ip = next_ip(ip)) {
      gen8_instruction *inst = &store[ip / 16];

      if (inst->opcode() == BRW_OPCODE_WHILE) {
         if (ip + inst->jip() <= start)
            return ip;
      }
   }
   assert(!"not reached");
   return start;
}

/* After program generation, go back and update the UIP and JIP of
 * BREAK, CONT, and HALT instructions to their correct locations.
 */
void
gen8_generator::patch_jump_targets()
{
   for (unsigned ip = 0; ip < next_inst_offset; ip = next_ip(ip)) {
      gen8_instruction *inst = &store[ip / 16];

      int block_end_ip = find_next_block_end(ip);
      switch (inst->opcode()) {
      case BRW_OPCODE_BREAK:
         assert(block_end_ip != 0);
         inst->set_jip(block_end_ip - ip);
         inst->set_uip(find_loop_end(ip) - ip);
         assert(inst->uip() != 0);
         assert(inst->jip() != 0);
         break;
      case BRW_OPCODE_CONTINUE:
         assert(block_end_ip != 0);
         inst->set_jip(block_end_ip - ip);
         inst->set_uip(find_loop_end(ip) - ip);
         assert(inst->uip() != 0);
         assert(inst->jip() != 0);
         break;
      case BRW_OPCODE_ENDIF:
         if (block_end_ip == 0)
            inst->set_jip(16);
         else
            inst->set_jip(block_end_ip - ip);
         break;
      case BRW_OPCODE_HALT:
         /* From the Sandy Bridge PRM (volume 4, part 2, section 8.3.19):
          *
          *    "In case of the halt instruction not inside any conditional
          *     code block, the value of <JIP> and <UIP> should be the
          *     same. In case of the halt instruction inside conditional code
          *     block, the <UIP> should be the end of the program, and the
          *     <JIP> should be end of the most inner conditional code block."
          *
          * The uip will have already been set by whoever set up the
          * instruction.
          */
         if (block_end_ip == 0) {
            inst->set_jip(inst->uip());
         } else {
            inst->set_jip(block_end_ip - ip);
         }
         assert(inst->uip() != 0);
         assert(inst->jip() != 0);
         break;
      }
   }
}

void
gen8_generator::DO()
{
   if (loop_stack_array_size < loop_stack_depth) {
      loop_stack_array_size *= 2;
      loop_stack = reralloc(mem_ctx, loop_stack, int, loop_stack_array_size);
   }
   loop_stack[loop_stack_depth++] = nr_inst;
}

gen8_instruction *
gen8_generator::BREAK()
{
   gen8_instruction *inst = next_inst(BRW_OPCODE_BREAK);
   inst->set_dst(retype(brw_null_reg(), BRW_REGISTER_TYPE_D));
   inst->set_src0(retype(brw_null_reg(), BRW_REGISTER_TYPE_D));
   inst->set_src1(brw_imm_d(0));
   inst->set_exec_size(default_state.exec_size);
   return inst;
}

gen8_instruction *
gen8_generator::CONTINUE()
{
   gen8_instruction *inst = next_inst(BRW_OPCODE_CONTINUE);
   inst->set_dst(brw_ip_reg());
   inst->set_src0(brw_ip_reg());
   inst->set_src1(brw_imm_d(0));
   inst->set_exec_size(default_state.exec_size);
   return inst;
}

gen8_instruction *
gen8_generator::WHILE()
{
   gen8_instruction *do_inst = &store[loop_stack[--loop_stack_depth]];
   gen8_instruction *while_inst = next_inst(BRW_OPCODE_WHILE);

   while_inst->set_dst(retype(brw_null_reg(), BRW_REGISTER_TYPE_D));
   while_inst->set_src0(retype(brw_null_reg(), BRW_REGISTER_TYPE_D));
   while_inst->set_src1(brw_imm_ud(0));
   while_inst->set_jip(16 * (do_inst - while_inst));
   while_inst->set_exec_size(default_state.exec_size);

   return while_inst;
}

gen8_instruction *
gen8_generator::HALT()
{
   gen8_instruction *inst = next_inst(BRW_OPCODE_HALT);
   inst->set_dst(retype(brw_null_reg(), BRW_REGISTER_TYPE_D));
   inst->set_src0(retype(brw_null_reg(), BRW_REGISTER_TYPE_D));
   inst->set_exec_size(default_state.exec_size);
   inst->set_mask_control(BRW_MASK_DISABLE);
   return inst;
}

void
gen8_generator::disassemble(FILE *out, int start, int end)
{
   bool dump_hex = false;

   for (int offset = start; offset < end; offset += 16) {
      gen8_instruction *inst = &store[offset / 16];
      printf("0x%08x: ", offset);

      if (dump_hex) {
         printf("0x%08x 0x%08x 0x%08x 0x%08x ",
                ((uint32_t *) inst)[3],
                ((uint32_t *) inst)[2],
                ((uint32_t *) inst)[1],
                ((uint32_t *) inst)[0]);
      }

      inst->disassemble(stdout, brw->gen);
   }
}
