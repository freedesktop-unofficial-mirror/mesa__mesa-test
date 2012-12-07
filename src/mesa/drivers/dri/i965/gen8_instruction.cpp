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

/** @file gen8_instruction.cpp
 *
 * A representation of a Gen8+ EU instruction, with helper methods to get
 * and set various fields.  This is the actual hardware format.
 */

#include "brw_defines.h"
#include "gen8_instruction.h"

void
gen8_instruction::set_dst(struct brw_reg reg)
{
   /* MRFs haven't existed since Gen7, so we better not be using them. */
   if (reg.file == BRW_MESSAGE_REGISTER_FILE) {
      reg.file = BRW_GENERAL_REGISTER_FILE;
      reg.nr += GEN7_MRF_HACK_START;
   }

   assert(reg.file != BRW_MESSAGE_REGISTER_FILE);

   if (reg.file == BRW_GENERAL_REGISTER_FILE)
      assert(reg.nr < BRW_MAX_GRF);

   set_dst_reg_file(reg.file);
   set_dst_reg_type(reg.type);

   assert(reg.address_mode == BRW_ADDRESS_DIRECT);

   set_dst_da_reg_nr(reg.nr);

   if (access_mode() == BRW_ALIGN_1) {
      /* Set Dst.SubRegNum[4:0] */
      set_dst_da1_subreg_nr(reg.subnr);

      /* Set Dst.HorzStride */
      if (reg.hstride == BRW_HORIZONTAL_STRIDE_0)
         reg.hstride = BRW_HORIZONTAL_STRIDE_1;
      set_dst_da1_hstride(reg.hstride);
   } else {
      /* Align16 SubRegNum only has a single bit (bit 4; bits 3:0 MBZ). */
      assert(reg.subnr == 0 || reg.subnr == 16);
      set_dst_da16_subreg_nr(reg.subnr >> 4);
      set_da16_writemask(reg.dw1.bits.writemask);
   }

   /* Generators should set a default exec_size of either 8 (SIMD4x2 or SIMD8)
    * or 16 (SIMD16), as that's normally correct.  However, when dealing with
    * small registers, we automatically reduce it to match the register size.
    */
   if (reg.width < BRW_EXECUTE_8)
      set_exec_size(reg.width);
}

void
gen8_instruction::validate_reg(struct brw_reg reg)
{
   int hstride_for_reg[] = {0, 1, 2, 4};
   int vstride_for_reg[] = {0, 1, 2, 4, 8, 16, 32, 64, 128, 256};
   int width_for_reg[] = {1, 2, 4, 8, 16};
   int execsize_for_reg[] = {1, 2, 4, 8, 16};
   int width, hstride, vstride, execsize;

   if (reg.file == BRW_IMMEDIATE_VALUE) {
      /* TODO(gen8): check immediate vectors */
      return;
   }

   if (reg.file == BRW_ARCHITECTURE_REGISTER_FILE)
      return;

   assert(reg.hstride >= 0 && reg.hstride < Elements(hstride_for_reg));
   hstride = hstride_for_reg[reg.hstride];

   if (reg.vstride == 0xf) {
      vstride = -1;
   } else {
      assert(reg.vstride >= 0 && reg.vstride < Elements(vstride_for_reg));
      vstride = vstride_for_reg[reg.vstride];
   }

   assert(reg.width >= 0 && reg.width < Elements(width_for_reg));
   width = width_for_reg[reg.width];

   assert(exec_size() >= 0 && exec_size() < Elements(execsize_for_reg));
   execsize = execsize_for_reg[exec_size()];

   /* Restrictions from 3.3.10: Register Region Restrictions. */
   /* 3. */
   assert(execsize >= width);

   /* 4. */
   if (execsize == width && hstride != 0) {
      assert(vstride == -1 || vstride == width * hstride);
   }

   /* 5. */
   if (execsize == width && hstride == 0) {
      /* no restriction on vstride. */
   }

   /* 6. */
   if (width == 1) {
      assert(hstride == 0);
   }

   /* 7. */
   if (execsize == 1 && width == 1) {
      assert(hstride == 0);
      assert(vstride == 0);
   }

   /* 8. */
   if (vstride == 0 && hstride == 0) {
      assert(width == 1);
   }

   /* 10. Check destination issues. */
}

void
gen8_instruction::set_src0(struct brw_reg reg)
{
   /* MRFs haven't existed since Gen7, so we better not be using them. */
   if (reg.file == BRW_MESSAGE_REGISTER_FILE) {
      reg.file = BRW_GENERAL_REGISTER_FILE;
      reg.nr += GEN7_MRF_HACK_START;
   }

   if (reg.file == BRW_GENERAL_REGISTER_FILE)
      assert(reg.nr < BRW_MAX_GRF);

   validate_reg(reg);

   set_src0_reg_file(reg.file);
   set_src0_reg_type(reg.type);
   set_src0_abs(reg.abs);
   set_src0_negate(reg.negate);

   assert(reg.address_mode == BRW_ADDRESS_DIRECT);

   if (reg.file == BRW_IMMEDIATE_VALUE) {
      data[3] = reg.dw1.ud;

      /* Required to set some fields in src1 as well: */
      set_src1_reg_file(0); /* arf */
      set_src1_reg_type(reg.type);
   } else {
      set_src0_da_reg_nr(reg.nr);

      if (access_mode() == BRW_ALIGN_1) {
         /* Set Src0.SubRegNum[4:0] */
         set_src0_da1_subreg_nr(reg.subnr);

         if (reg.width == BRW_WIDTH_1 && exec_size() == BRW_EXECUTE_1) {
            set_src0_da1_hstride(BRW_HORIZONTAL_STRIDE_0);
            set_src0_vert_stride(BRW_VERTICAL_STRIDE_0);
         } else {
            set_src0_da1_hstride(reg.hstride);
            set_src0_vert_stride(reg.vstride);
         }
         set_src0_da1_width(reg.width);

      } else {
         /* Align16 SubRegNum only has a single bit (bit 4; bits 3:0 MBZ). */
         assert(reg.subnr == 0 || reg.subnr == 16);
         set_src0_da16_subreg_nr(reg.subnr >> 4);

         set_src0_da16_swiz_x(BRW_GET_SWZ(reg.dw1.bits.swizzle, BRW_CHANNEL_X));
         set_src0_da16_swiz_y(BRW_GET_SWZ(reg.dw1.bits.swizzle, BRW_CHANNEL_Y));
         set_src0_da16_swiz_z(BRW_GET_SWZ(reg.dw1.bits.swizzle, BRW_CHANNEL_Z));
         set_src0_da16_swiz_w(BRW_GET_SWZ(reg.dw1.bits.swizzle, BRW_CHANNEL_W));

         /* This is an oddity of the fact that we're using the same
          * descriptions for registers in both Align16 and Align1 modes.
          */
         if (reg.vstride == BRW_VERTICAL_STRIDE_8)
            set_src0_vert_stride(BRW_VERTICAL_STRIDE_4);
         else
            set_src0_vert_stride(reg.vstride);
      }
   }
}

void
gen8_instruction::set_src1(struct brw_reg reg)
{
   /* MRFs haven't existed since Gen7, so we better not be using them. */
   if (reg.file == BRW_MESSAGE_REGISTER_FILE) {
      reg.file = BRW_GENERAL_REGISTER_FILE;
      reg.nr += GEN7_MRF_HACK_START;
   }

   if (reg.file == BRW_GENERAL_REGISTER_FILE)
      assert(reg.nr < BRW_MAX_GRF);

   validate_reg(reg);

   set_src1_reg_file(reg.file);
   set_src1_reg_type(reg.type);
   set_src1_abs(reg.abs);
   set_src1_negate(reg.negate);

   /* Only src1 can be an immediate in two-argument instructions. */
   assert(src0_reg_file() != BRW_IMMEDIATE_VALUE);

   assert(reg.address_mode == BRW_ADDRESS_DIRECT);

   if (reg.file == BRW_IMMEDIATE_VALUE) {
      data[3] = reg.dw1.ud;
   } else {
      set_src1_da_reg_nr(reg.nr);

      if (access_mode() == BRW_ALIGN_1) {
         /* Set Src0.SubRegNum[4:0] */
         set_src1_da1_subreg_nr(reg.subnr);

         if (reg.width == BRW_WIDTH_1 && exec_size() == BRW_EXECUTE_1) {
            set_src1_da1_hstride(BRW_HORIZONTAL_STRIDE_0);
            set_src1_vert_stride(BRW_VERTICAL_STRIDE_0);
         } else {
            set_src1_da1_hstride(reg.hstride);
            set_src1_vert_stride(reg.vstride);
         }
         set_src1_da1_width(reg.width);
      } else {
         /* Align16 SubRegNum only has a single bit (bit 4; bits 3:0 MBZ). */
         assert(reg.subnr == 0 || reg.subnr == 16);
         set_src1_da16_subreg_nr(reg.subnr >> 4);

         set_src1_da16_swiz_x(BRW_GET_SWZ(reg.dw1.bits.swizzle, BRW_CHANNEL_X));
         set_src1_da16_swiz_y(BRW_GET_SWZ(reg.dw1.bits.swizzle, BRW_CHANNEL_Y));
         set_src1_da16_swiz_z(BRW_GET_SWZ(reg.dw1.bits.swizzle, BRW_CHANNEL_Z));
         set_src1_da16_swiz_w(BRW_GET_SWZ(reg.dw1.bits.swizzle, BRW_CHANNEL_W));

         /* This is an oddity of the fact that we're using the same
          * descriptions for registers in both Align16 and Align1 modes.
          */
         if (reg.vstride == BRW_VERTICAL_STRIDE_8)
            set_src1_vert_stride(BRW_VERTICAL_STRIDE_4);
         else
            set_src1_vert_stride(reg.vstride);
      }
   }
}

/**
 * Set the Message Descriptor and Extended Message Descriptor fields
 * for SEND messages.
 *
 * \note This zeroes out the Function Control bits, so it must be called
 *       \b before filling out any message-specific data.  Callers can
 *       choose not to fill in irrelevant bits; they will be zero.
 */
void
gen8_instruction::set_message_descriptor(enum brw_message_target sfid,
                                         unsigned msg_length,
                                         unsigned response_length,
                                         bool header_present,
                                         bool end_of_thread)
{
   set_src1(brw_imm_d(0));

   set_sfid(sfid);
   set_mlen(msg_length);
   set_rlen(response_length);
   set_header_present(header_present);
   set_eot(end_of_thread);
}

void
gen8_instruction::set_urb_message(enum brw_urb_write_flags flags,
                                  unsigned msg_length,
                                  unsigned response_length,
                                  unsigned offset,
                                  bool interleave)
{
   set_message_descriptor(BRW_SFID_URB, msg_length, response_length,
                          true, flags & BRW_URB_WRITE_EOT);
   set_src0(brw_vec8_grf(GEN7_MRF_HACK_START + 1, 0));
   if (flags & BRW_URB_WRITE_OWORD) {
      assert(msg_length == 2);
      set_urb_opcode(BRW_URB_OPCODE_WRITE_OWORD);
   } else {
      set_urb_opcode(BRW_URB_OPCODE_WRITE_HWORD);
   }
   set_urb_global_offset(offset);
   set_urb_interleave(interleave);
   set_urb_per_slot_offset(flags & BRW_URB_WRITE_PER_SLOT_OFFSET ? 1 : 0);
}

void
gen8_instruction::set_sampler_message(unsigned binding_table_index,
                                      unsigned sampler,
                                      unsigned msg_type,
                                      unsigned response_length,
                                      unsigned msg_length,
                                      bool header_present,
                                      unsigned simd_mode)
{
   set_message_descriptor(BRW_SFID_SAMPLER, msg_length, response_length,
                          header_present, false);

   set_binding_table_index(binding_table_index);
   set_sampler(sampler);
   set_sampler_msg_type(msg_type);
   set_sampler_simd_mode(simd_mode);
}

void
gen8_instruction::set_dp_message(enum brw_message_target sfid,
                                 unsigned binding_table_index,
                                 unsigned msg_type,
                                 unsigned msg_control,
                                 unsigned mlen,
                                 unsigned rlen,
                                 bool header_present,
                                 bool end_of_thread)
{
   /* Binding table index is from 0..255 */
   assert((binding_table_index & 0xff) == binding_table_index);

   /* Message Type is only 5 bits */
   assert((msg_type & 0x1f) == msg_type);

   /* Message Control is only 6 bits */
   assert((msg_control & 0x3f) == msg_control);

   set_message_descriptor(sfid, mlen, rlen, header_present, end_of_thread);
   set_function_control(binding_table_index | msg_type << 14 | msg_control << 8);
}
