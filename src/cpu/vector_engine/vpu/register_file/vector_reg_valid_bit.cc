/*
 * Copyright (c) 2020 Barcelona Supercomputing Center
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met: redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer;
 * redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution;
 * neither the name of the copyright holders nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Author: Cristóbal Ramírez
 */

#include "cpu/vector_engine/vpu/register_file/vector_reg_valid_bit.hh"

#include "debug/VectorValidBit.hh"

/**
 * Valid bits
 */
VectorValidBit::VectorValidBit(VectorValidBitParams *p):
SimObject(p), PhysicalRegs(p->PhysicalRegs), PhyReg_scalar(32),
AdditionalRegs(5)
{
    for (uint64_t i = 0; i < 32; i++) {
        reg_valid_bit.push_back(1);
    }
    for (uint64_t i = 32; i < PhysicalRegs; i++) {
        reg_valid_bit.push_back(0);
        //vector freelist only contains 32~PhysicalRegs
    }
    for (uint64_t i = 0; i < PhyReg_scalar; i++) {
        scalar_reg_valid_bit.push_back(1);
    }
    for (uint64_t i = PhyReg_scalar;
        i < PhyReg_scalar + AdditionalRegs; i++) {
        scalar_reg_valid_bit.push_back(0);
    }
}

VectorValidBit::~VectorValidBit()
{
}

void
VectorValidBit::set_preg_valid_bit(int idx, int val)
{
    assert(idx <= PhysicalRegs);
    reg_valid_bit[idx] = val;
    DPRINTF(VectorValidBit, "Setting the validbit reg %d with %d\n"
            ,idx,val);
}

int
VectorValidBit::get_preg_valid_bit(int idx)
{
    assert((idx <= PhysicalRegs));
    return reg_valid_bit[idx];
}

void
VectorValidBit::set_pscalar_reg_valid_bit(int idx, int val)
{
    assert(idx <= PhysicalRegs);

    reg_valid_bit[idx] = val;
    DPRINTF(VectorValidBit, "Setting the validbit reg %d with %d\n"
            ,idx,val);
}

int
VectorValidBit::get_pscalar_reg_valid_bit(int idx)
{
    assert((idx <= PhysicalRegs));
    return reg_valid_bit[idx];
}

void
VectorValidBit::print_valid_bit()
{

}

VectorValidBit *
VectorValidBitParams::create()
{
    return new VectorValidBit(this);
}