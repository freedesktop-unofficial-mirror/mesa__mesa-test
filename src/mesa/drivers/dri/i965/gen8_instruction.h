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

/** @file gen8_instruction.h
 *
 * A representation of a Gen8+ EU instruction, with helper methods to get
 * and set various fields.  This is the actual hardware format.
 */

#pragma once

extern "C" {
#include "main/macros.h"
} /* extern "C" */

#include "brw_reg.h"

struct gen8_instruction {
   uint32_t data[4];

#define F(name, high, low) \
   inline void set_##name(unsigned v) { set_bits(high, low, v); } \
   inline unsigned name() const { return bits(high, low); }

   /**
    * Direct addressing only:
    *  @{
    */
   F(src1_da_reg_nr,      108, 101);
   F(src0_da_reg_nr,       76,  69);
   F(dst_da1_hstride,      62,  61);
   F(dst_da_reg_nr,        60,  53);
   F(dst_da16_subreg_nr,   52,  52);
   F(dst_da1_subreg_nr,    52,  48);
   F(da16_writemask,       51,  48); /* Dst.ChanEn */
   /** @} */

   F(src1_vert_stride,    120, 117)
   F(src1_da1_width,      116, 114)
   F(src1_da16_swiz_w,    115, 114)
   F(src1_da16_swiz_z,    113, 112)
   F(src1_da1_hstride,    113, 112)
   F(src1_address_mode,   111, 111)
   /** Src1.SrcMod @{ */
   F(src1_negate,         110, 110)
   F(src1_abs,            109, 109)
   /** @} */
   F(src1_da16_subreg_nr, 100, 100)
   F(src1_da1_subreg_nr,  100,  96)
   F(src1_da16_swiz_y,     99,  98)
   F(src1_da16_swiz_x,     97,  96)
   F(src1_reg_type,        94,  91)
   F(src1_reg_file,        90,  89)
   F(src0_vert_stride,     88,  85)
   F(src0_da1_width,       84,  82)
   F(src0_da16_swiz_w,     83,  82)
   F(src0_da16_swiz_z,     81,  80)
   F(src0_da1_hstride,     81,  80)
   F(src0_address_mode,    79,  79)
   /** Src0.SrcMod @{ */
   F(src0_negate,          78,  78)
   F(src0_abs,             77,  77)
   /** @} */
   F(src0_da16_subreg_nr,  68,  68)
   F(src0_da1_subreg_nr,   68,  64)
   F(src0_da16_swiz_y,     67,  66)
   F(src0_da16_swiz_x,     65,  64)
   F(dst_address_mode,     63,  63)
   F(src0_reg_type,        46,  43)
   F(src0_reg_file,        42,  41)
   F(dst_reg_type,         40,  37)
   F(dst_reg_file,         36,  35)
   F(mask_control,         34,  34)
   F(flag_reg_nr,          33,  33)
   F(flag_subreg_nr,       32,  32)
   F(saturate,             31,  31)
   F(branch_control,       30,  30)
   F(debug_control,        30,  30)
   F(cmpt_control,         29,  29)
   F(acc_wr_control,       28,  28)
   F(cond_modifier,        27,  24)
   F(exec_size,            23,  21)
   F(pred_inv,             20,  20)
   F(pred_control,         19,  16)
   F(thread_control,       15,  14)
   F(qtr_control,          13,  12)
   F(nib_control,          11,  11)
   F(no_dd_check,          10,  10)
   F(no_dd_clear,           9,   9)
   F(access_mode,           8,   8)
   /* Bit 7 is Reserve d (for future Opcode expansion) */
   F(opcode,                6,   0)

   /**
    * Three-source instructions:
    *  @{
    */
   F(src2_3src_reg_nr,    125, 118)
   F(src2_3src_subreg_nr, 117, 115)
   F(src2_3src_swizzle,   114, 107)
   F(src2_3src_rep_ctrl,  106, 106)
   F(src1_3src_reg_nr,    104,  97)
   F(src1_3src_subreg_hi,  96,  96)
   F(src1_3src_subreg_lo,  95,  94)
   F(src1_3src_swizzle,    93,  86)
   F(src1_3src_rep_ctrl,   85,  85)
   F(src0_3src_reg_nr,     83,  76)
   F(src0_3src_subreg_nr,  75,  73)
   F(src0_3src_swizzle,    72,  65)
   F(src0_3src_rep_ctrl,   64,  64)
   F(dst_3src_reg_nr,      63,  56)
   F(dst_3src_subreg_nr,   55,  53)
   F(dst_3src_writemask,   52,  49)
   F(dst_3src_type,        48,  46)
   F(src_3src_type,        45,  43)
   F(src2_3src_negate,     42,  42)
   F(src2_3src_abs,        41,  41)
   F(src1_3src_negate,     40,  40)
   F(src1_3src_abs,        39,  39)
   F(src0_3src_negate,     38,  38)
   F(src0_3src_abs,        37,  37)
   /** @} */

   /**
    * Fields for SEND messages:
    *  @{
    */
   F(eot,                 127, 127)
   F(mlen,                124, 121)
   F(rlen,                120, 116)
   F(header_present,      115, 115)
   F(function_control,    114,  96)
   F(sfid,                 27,  24)
   F(math_function,        27,  24)
   /** @} */

   /**
    * URB message function control bits:
    *  @{
    */
   F(urb_per_slot_offset, 113, 113)
   F(urb_interleave,      111, 111)
   F(urb_global_offset,   110, 100)
   F(urb_opcode,           99,  96)
   /** @} */

   /**
    * Sampler message function control bits:
    *  @{
    */
   F(sampler_simd_mode,   114, 113)
   F(sampler_msg_type,    112, 108)
   F(sampler,             107, 104)
   F(binding_table_index, 103,  96)
   /** @} */
#undef F

   /**
    * Flow control instruction bits:
    *  @{
    */
   inline unsigned uip() const { return data[2]; }
   inline void set_uip(unsigned uip) { data[2] = uip; }
   inline unsigned jip() const { return data[3]; }
   inline void set_jip(unsigned jip) { data[3] = jip; }
   /** @} */

   inline int src1_imm_d() const { return data[3]; }
   inline unsigned src1_imm_ud() const { return data[3]; }
   inline float src1_imm_f() const { fi_type ft; ft.u = data[3]; return ft.f; }

   void set_dst(struct brw_reg reg);
   void set_src0(struct brw_reg reg);
   void set_src1(struct brw_reg reg);

   void set_urb_message(enum brw_urb_write_flags, unsigned mlen, unsigned rlen,
                        unsigned offset, bool interleave);

   void set_sampler_message(unsigned binding_table_index, unsigned sampler,
                            unsigned msg_type, unsigned rlen, unsigned mlen,
                            bool header_present, unsigned simd_mode);

   void set_dp_message(enum brw_message_target sfid,
                       unsigned binding_table_index,
                       unsigned msg_type,
                       unsigned msg_control,
                       unsigned msg_length,
                       unsigned response_length,
                       bool header_present,
                       bool end_of_thread);

private:
   inline unsigned bits(unsigned high, unsigned low) const;
   inline void set_bits(unsigned high, unsigned low, unsigned value);

   void validate_reg(struct brw_reg reg);

   void set_message_descriptor(enum brw_message_target sfid,
                               unsigned msg_length,
                               unsigned response_length,
                               bool header_present,
                               bool end_of_thread);
};

/**
 * Fetch a set of contiguous bits from the instruction.
 *
 * Bits indexes range from 0..127; fields may not cross 32-bit boundaries.
 */
inline unsigned
gen8_instruction::bits(unsigned high, unsigned low) const
{
   /* We assume the field doesn't cross 32-bit boundaries. */
   const unsigned word = high / 32;
   assert(word == low / 32);

   high %= 32;
   low %= 32;

   const unsigned mask = (((1 << (high - low + 1)) - 1) << low);

   return (data[word] & mask) >> low;
}

/**
 * Set bits in the instruction, with proper shifting and masking.
 *
 * Bits indexes range from 0..127; fields may not cross 32-bit boundaries.
 */
inline void
gen8_instruction::set_bits(unsigned high, unsigned low, unsigned value)
{
   const unsigned word = high / 32;
   assert(word == low / 32);

   high %= 32;
   low %= 32;

   const unsigned mask = (((1 << (high - low + 1)) - 1) << low);

   data[word] = (data[word] & ~mask) | ((value << low) & mask);
}
