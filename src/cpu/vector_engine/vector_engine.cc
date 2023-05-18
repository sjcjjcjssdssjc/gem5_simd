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

#include "cpu/vector_engine/vector_engine.hh"

#include <cassert>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#include "base/logging.hh"
#include "base/types.hh"
#include "cpu/translation.hh"
#include "debug/VectorEngine.hh"
#include "debug/VectorEngineInfo.hh"
#include "debug/VectorInst.hh"
#include "sim/faults.hh"
#include "sim/sim_object.hh"

VectorEngine::VectorEngine(VectorEngineParams *p) :
SimObject(p),
vector_config(p->vector_config),
VectorCacheMasterId(p->system->getMasterId(this, name() + ".vector_cache")),
vectormem_port(name() + ".vector_mem_port", this, p->vector_rf_ports),
vector_reg(p->vector_reg),
uniqueReqId(0),
num_clusters(p->num_clusters),
num_lanes(p->num_lanes),
vector_rob(p->vector_rob),
vector_lane(p->vector_lane),
vector_memory_unit(p->vector_memory_unit),
vector_inst_queue(p->vector_inst_queue),
vector_rename(p->vector_rename),
vector_reg_validbit(p->vector_reg_validbit),
last_vtype(0),
last_vl(0)
{
    //create independent ports
    for (uint8_t i=0; i < p->vector_rf_ports; ++i) {
        VectorRegMasterIds.push_back(p->system->getMasterId(this, name()
            + ".vector_reg" + std::to_string(i)));
        VectorRegPorts.push_back(VectorRegPort(name()
            + ".vector_reg_port", this, i));
    }

    DPRINTF(VectorEngineInfo,"VPU created:\n");
    DPRINTF(VectorEngineInfo,"Vector Renaming Enabled\n");
    DPRINTF(VectorEngineInfo,"Number of Physical Registers: %d \n"
        ,vector_rename->PhysicalRegs );
    DPRINTF(VectorEngineInfo,"Maximum VL: %d-bits\n"
        , vector_config->get_mvl_lmul1_bits() );
    DPRINTF(VectorEngineInfo,"Register File size: %dKB\n"
        , (float)vector_reg->get_size()/1024.0 );
    DPRINTF(VectorEngineInfo,"Vector Reorder buffer size: %d\n"
        , vector_rob->ROB_Size );
    DPRINTF(VectorEngineInfo,"Vector Arithmetic Queue Size: %d \n"
        ,vector_inst_queue->vector_arith_queue_size );
    DPRINTF(VectorEngineInfo,"Vector Memory Queue Size: %d \n"
        ,vector_inst_queue->vector_mem_queue_size );
    if (vector_inst_queue->OoO_queues)
        DPRINTF(VectorEngineInfo,"Out-of-Order Issue Logic\n");
    else
        DPRINTF(VectorEngineInfo,"In-Order Issue Logic\n");
    if (num_lanes<num_clusters) {
        panic("Invalid VPU Configuration\n");
    }
    DPRINTF(VectorEngineInfo,"Number of Clusters: %d\n",num_clusters);
    DPRINTF(VectorEngineInfo,"Lanes per Cluster: %d\n",num_lanes/num_clusters);
    DPRINTF(VectorEngineInfo,"Vector Memory Unit Port connected to l2 bus\n");
}

VectorEngine::~VectorEngine()
{
}


void
VectorEngine::regStats()
{
    SimObject::regStats();

    VectorArithmeticIns
        .name(name() + ".VectorArithmeticIns")
        .desc("Number vector arithmetic instructions");
    VectorMemIns
        .name(name() + ".VectorMemIns")
        .desc("Number vector arithmetic instructions");
    VectorConfigIns
        .name(name() + ".VectorConfigIns")
        .desc("Number vector configuration instructions");
    SumVL
        .name(name() + ".SumVL")
        .desc("Sum of all instruction's vector length to obtain the average ");
    AverageVL
        .name(name() + ".AverageVL")
        .desc("SumVL divided by the number of instructions")
        .precision(1);
    AverageVL = SumVL/(VectorMemIns+VectorArithmeticIns);

}


bool
VectorEngine::requestGrant(RiscvISA::VectorStaticInst *insn)
{
    bool rob_entry_available = !vector_rob->rob_full();
    bool queue_slots_available = ((insn->isVectorInstArith()
        || insn->isVecConfig()) && !vector_inst_queue->arith_queue_full())
        || (insn->isVectorInstMem()
        && !vector_inst_queue->mem_queue_full());

    /*
     * Usually, the Vector engine must ensure to have at least 1 physical register to accept an instruction.
     * However, for widening operations this is different. When Widening operation with LMUL=2 means that
     * the vector engine must have at least 2 physical registers available in order to acceept the instruction.
     * When LMUL=4, means that the vector engine must have at least 4 physical registers available in order 
     * to accept the instruction. And finally when LMUL=8, means that the vector engine must have at least 
     * 8 physical registers available in order to accept the instruction.
     */ 
    bool mask_dst = insn->isMaskDst();

    bool enough_physical_regs = ((last_lmul == 1) || mask_dst) ? vector_rename->frl_elements() >= 1 :
                                (last_lmul == 2) ? vector_rename->frl_elements() >= 2 :
                                (last_lmul == 4) ? vector_rename->frl_elements() >= 4 :
                                (last_lmul == 8) ? vector_rename->frl_elements() >= 8 : 0;
    //DPRINTF(VectorInst,"rob_entry_available %d, queue_slots_available %d, enough_physical_regs %d\n",
    //                            rob_entry_available,queue_slots_available,enough_physical_regs);
    return  enough_physical_regs && queue_slots_available
        && rob_entry_available;
    //return  !vector_rename->frl_empty() && queue_slots_available
    //    && rob_entry_available;
}

bool
VectorEngine::isOccupied()
{
    bool cluster_ocuppied = 0;
    for (int i=0 ; i< num_clusters ; i++)
    {
        cluster_ocuppied = cluster_ocuppied || vector_lane[i]->isOccupied();
    }

    bool data_in_queues = 0;
    data_in_queues = (vector_inst_queue->Instruction_Queue.size() > 0)
        || (vector_inst_queue->Memory_Queue.size() > 0);

    bool rob_empty = vector_rob->rob_empty();

    return (cluster_ocuppied || data_in_queues
        || vector_memory_unit->isOccupied()
        || vector_inst_queue->isOccupied() || !rob_empty);
}

bool
VectorEngine::cluster_available()
{
    bool available = 0;

    for (int i = 0; i < num_clusters; i++)
    {
        available = available || !vector_lane[i]->isOccupied();
    }

    return available;
}

void
VectorEngine::printConfigInst(RiscvISA::VectorStaticInst& insn, uint64_t src1,uint64_t src2)
{
    uint64_t pc = insn.getPC();
    DPRINTF(VectorInst,"inst: %s vl:%d, vtype:%d PC 0x%X\n"
        ,insn.getName(),src1,src2,*(uint64_t*)&pc );
}

void
VectorEngine::printMemInst(RiscvISA::VectorStaticInst& insn,VectorDynInst *vector_dyn_insn)
{

}

void
VectorEngine::printArithInst(RiscvISA::VectorStaticInst& insn,VectorDynInst *vector_dyn_insn)
{

}

void
VectorEngine::renameVectorInst(RiscvISA::VectorStaticInst& insn, VectorDynInst *vector_dyn_insn)
{
    /*
     * The bigger LMUL value is 8, then for LMUL=8 it is needed to assign 8 physical registers
     * for LMUL = 4 it is needed to assign 4 physical registers
     * for LMUL = 2 it is needed to assign 2 physical registers
     * for LMUL = 1 it is needed to assign 1 physical register
     */
    uint64_t vd;
    uint64_t vs1,vs2,vs3;
    uint64_t rs1,rs2;
    uint64_t Prs1,Prs2;

    vd = insn.vd();
    vs1 = insn.vs1();
    vs2 = insn.vs2();
    vs3 = insn.vs3();

    rs1 = insn.rs1();
    rs2 = insn.rs2();

    masked_op = (insn.vm() == 0);

    vx_op = (insn.func3() == 4) || (insn.func3() == 6);
    vf_op = (insn.func3() == 5);
    vi_op = (insn.func3() == 3);

    uint8_t mop = insn.mop();
    bool indexed = (mop == 3);
    bool strided = (mop == 2);

    if (insn.isVectorInstMem()) {
        // TD: masked and strided memory operations are not implemented
        if (insn.isLoad()) {
            Prs1 = Prs2 = 1024;
            if (strided) {
                Prs1 = vector_rename->get_preg_ratscalar(rs1);
                Prs2 = vector_rename->get_preg_ratscalar(rs2);
                DPRINTF(VectorRename, "strided load Prs1: %d Prs2: %d\n",
                Prs1, Prs2);
            }
            PMask = masked_op ? vector_rename->get_preg_rat(0) : 1024;
            Pvs2  = indexed ? vector_rename->get_preg_rat(vs2) : 1024;
            PDst  = vector_rename->get_frl();
            POldDst = vector_rename->get_preg_rat(vd);

            vector_dyn_insn->set_renamed_src1(strided ? Prs1 : 1024);
            vector_dyn_insn->set_renamed_src2(indexed ? Pvs2 : Prs2);

            vector_dyn_insn->set_renamed_mask(PMask);
            vector_dyn_insn->set_renamed_dst(PDst);
            vector_dyn_insn->set_renamed_old_dst(POldDst);

            vector_rename->set_preg_rat(vd, PDst);
            vector_reg_validbit->set_preg_valid_bit(PDst, 0);
        } else if (insn.isStore()) {
            Prs1 = Prs2 = 1024;
            if (strided) {
                Prs1 = vector_rename->get_preg_ratscalar(rs1);
                Prs2 = vector_rename->get_preg_ratscalar(rs2);
                DPRINTF(VectorRename, "strided store Prs1: %d Prs2: %d\n",
                Prs1, Prs2);
            }
            Pvs2 = (indexed) ? vector_rename->get_preg_rat(vs2) : 1024;
            Pvs3 = vector_rename->get_preg_rat(vs3);
            PMask = masked_op ? vector_rename->get_preg_rat(0) : 1024;
            if (strided) {
                vector_dyn_insn->set_renamed_src1(Prs1);
                vector_dyn_insn->set_renamed_src2(Prs2);
            } else {
                vector_dyn_insn->set_renamed_src2(Pvs2);
            }
            /* stored data */
            vector_dyn_insn->set_renamed_src3(Pvs3);
            vector_dyn_insn->set_renamed_mask(PMask);
        }
    }
    else if (insn.isVectorInstArith()) {
        PMask = masked_op ? vector_rename->get_preg_rat(0) : 1024;
        PDst = (dst_write_ena) ? vector_rename->get_frl() : 1024;
        POldDst = (dst_write_ena) ? vector_rename->get_preg_rat(vd) : 1024;
        Pvs2 = vector_rename->get_preg_rat(vs2);

        vector_dyn_insn->set_renamed_mask(PMask);
        vector_dyn_insn->set_renamed_dst(PDst);
        vector_dyn_insn->set_renamed_old_dst(POldDst);
        vector_dyn_insn->set_renamed_src2(Pvs2);

        if (!(vx_op || vf_op || vi_op) && !insn.arith1Src()) {
            // Physical source 1
            Pvs1 = vector_rename->get_preg_rat(vs1);
            vector_dyn_insn->set_renamed_src1(Pvs1);
        } else if (vx_op || vf_op || vi_op) {
            Prs1 = vector_rename->get_preg_ratscalar(rs1);
            vector_dyn_insn->set_renamed_src1(Prs1);
        }

        // dst_write_ena set when instruction has a vector destination register
        if (dst_write_ena) {
            // Setting the New Destination in the RAT structure
            // Setting to 0 the new physical destinatio valid bit
            vector_rename->set_preg_rat(vd,PDst);
            vector_reg_validbit->set_preg_valid_bit(PDst, 0);
        }
    } else {
        panic("Invalid Vector Instruction insn=%#h\n", insn.machInst);
    }
}

void
VectorEngine::dispatch(RiscvISA::VectorStaticInst& insn, ExecContextPtr& xc,
    uint64_t src1,uint64_t src2,std::function<void()> dependencie_callback) {
    if (insn.isVecConfig()) {
        last_vtype = src2;
        last_vl = src1;
        dependencie_callback();
        printConfigInst(insn,src1,src2);
        VectorConfigIns++;
        DPRINTF(VectorEngine,"Settig vl %d , sew %d , lmul %d\n",last_vl,vector_config->get_vtype_sew(last_vtype),vector_config->get_vtype_lmul(last_vtype));
        return;
    }

    // There is a very initial support for widening convert from int32 to fp64
    // Basic support for all widening and narrowing conversions will be addded,
    // but limiting to MVL/2 in case of widining, ince we are not implementing register grouping.
    // This will allow to convert between int32,int64,fp32 and fp64.
    if(/*insn.isWidening() ||*/ insn.isNarrowing()){
        panic("Widening/Narrowing vector instructions are not fully suported \n");
    }

    last_lmul = vector_config->get_vtype_lmul(last_vtype);
    if (last_lmul > 1) {
        panic("LMUL>1 is not suported \n");
    }

    //DPRINTF(VectorInst,"src1 %d, src2 %d\n",src1,src2);

    /* Be sure that the instruction was added to some group in base.isa */
    if (insn.isVectorInstArith()) {
        assert(insn.arith1Src() | insn.arith2Srcs() | insn.arith3Srcs());
    }

    if ((vector_inst_queue->Instruction_Queue.size()==0)
        && (vector_inst_queue->Memory_Queue.size()==0)) {
        vector_inst_queue->startTicking(*this/*,dependencie_callback*/);
    }

    dst_write_ena = !insn.VectorToScalar();

    VectorDynInst *vector_dyn_insn = new VectorDynInst();

    renameVectorInst(insn,vector_dyn_insn);

    if (vector_rob->rob_empty()) {
        vector_rob->startTicking(*this);
    }

    if (insn.isVectorInstMem()) {
        dependencie_callback();
        uint32_t rob_entry = vector_rob->
        set_rob_entry(vector_dyn_insn->get_renamed_old_dst(),
                      insn.isLoad(), 0);
        vector_dyn_insn->set_rob_num(rob_entry);
        vector_inst_queue->Memory_Queue.push_back(
            new InstQueue::QueueEntry(insn,vector_dyn_insn,xc,
                NULL,src1,src2,last_vtype,last_vl));
        printMemInst(insn,vector_dyn_insn);
    }
    else if (insn.isVectorInstArith()) {
        if (dst_write_ena) {
            dependencie_callback();
            bool to_scalar = insn.VectorToScalar();
            uint32_t rob_entry = vector_rob->set_rob_entry(
                vector_dyn_insn->get_renamed_old_dst(), 1,
                to_scalar);
            vector_dyn_insn->set_rob_num(rob_entry);
            vector_inst_queue->Instruction_Queue.push_back(
            new InstQueue::QueueEntry(insn,vector_dyn_insn,xc,
            NULL,src1,src2,last_vtype,last_vl));
        } else {
            uint32_t rob_entry = vector_rob->set_rob_entry(0, 0, 0);
            vector_dyn_insn->set_rob_num(rob_entry);
            vector_inst_queue->Instruction_Queue.push_back(
            new InstQueue::QueueEntry(insn,vector_dyn_insn,xc,
            dependencie_callback,src1,src2,last_vtype,last_vl));
        }
        printArithInst(insn,vector_dyn_insn);
    } else {
        panic("Invalid Vector Instruction, insn=%X\n", insn.machInst);
    }
}

void
VectorEngine::issue(RiscvISA::VectorStaticInst& insn,VectorDynInst *dyn_insn,
    ExecContextPtr& xc,uint64_t src1,uint64_t src2,uint64_t vtype,
    uint64_t vl, std::function<void(Fault fault)> done_callback) {
    uint64_t pc = insn.getPC();
    
    if (insn.isVectorInstMem()) {
        // There is no scalar in this
        VectorMemIns++;
        DPRINTF(VectorEngine, "Sending instruction %s to VMU, pc 0x%lx\n"
            ,insn.getName(), *(uint64_t*)& pc);
        // done_callback is write back.
        vector_memory_unit->issue(*this,insn,dyn_insn,xc,src1,src2,vtype,
            vl, done_callback);
        SumVL = SumVL.value() + vector_config->vector_length_in_bits(vl,vtype);
    }
    else if (insn.isVectorInstArith()) {
        SumVL = SumVL.value() + vector_config->vector_length_in_bits(vl,vtype);
        VectorArithmeticIns++;
        uint8_t lane_id_available = 0;
        for (int i = 0 ; i < num_clusters; i++) {
            if (!vector_lane[i] -> isOccupied()) {
                lane_id_available = i;
            }
        }
        DPRINTF(VectorEngine, "Sending instruction %s to Lanes, pc 0x%lx\n",
            insn.getName(), *(uint64_t*)&pc);
        vector_lane[lane_id_available]->issue(*this,insn,dyn_insn,xc,src1,
            vtype, vl,done_callback);
    } else {
        panic("Invalid Vector Instruction, insn=%#h\n", insn.machInst);
    }
}

Port &
VectorEngine::getPort(const std::string &if_name, PortID idx)
{
    if (if_name == "vector_reg_port")
        return VectorRegPorts[idx];
    else if (if_name == "vector_mem_port")
        return getVectorMemPort();
    else
        return getPort(if_name, idx);
}

VectorEngine::VectorMemPort::VectorMemPort(const std::string& name,
    VectorEngine* owner, uint8_t channels) :
    MasterPort(name, owner), owner(owner)
{
    //create the queues for each of the channels to the Vector Cache
    for (uint8_t i = 0; i<channels; ++i) {
        laCachePktQs.push_back(std::deque<PacketPtr>());
    }
}

VectorEngine::VectorMemPort::~VectorMemPort()
{
}

VectorEngine::VectorMemPort::Tlb_Translation::Tlb_Translation(
    VectorEngine *owner):
    event(this, true), owner(owner)
{
}

VectorEngine::VectorMemPort::Tlb_Translation::~Tlb_Translation()
{
}

void
VectorEngine::VectorMemPort::Tlb_Translation::finish(const Fault &_fault,
    const RequestPtr &_req, ThreadContext *_tc, BaseTLB::Mode _mode)
{
    fault = _fault;
}

void
VectorEngine::VectorMemPort::Tlb_Translation::finish(const Fault _fault,
    uint64_t latency)
{
    fault = _fault;
}

void
VectorEngine::VectorMemPort::Tlb_Translation::markDelayed()
{
    panic("Tlb_Translation::markDelayed not implemented");
}

std::string
VectorEngine::VectorMemPort::Tlb_Translation::name()
{
    //TODO: IDK what to put here
    return "Tlb_Translation";
}

void
VectorEngine::VectorMemPort::Tlb_Translation::translated()
{
}

// submits a request for translation to the execCacheTlb
// MEMCPY data if present, so caller must deallocate it himself!
bool
VectorEngine::VectorMemPort::startTranslation(Addr addr, uint8_t *data,
    uint64_t size, BaseTLB::Mode mode, ThreadContext *tc, uint64_t req_id,
    uint8_t channel)
{
    Process * p = tc->getProcessPtr();
    Addr page1 = p->pTable->pageAlign(addr);
    Addr page2 = p->pTable->pageAlign(addr+size-1);
    assert(page1 == page2);

    //NOTE: need to make a buffer for reads so cache can write to it!
    uint8_t *ndata = new uint8_t[size];
    if (data != nullptr) {
        assert(mode == BaseTLB::Write);
        memcpy(ndata, data, size);
    } else {
      //put a pattern here for debugging
      memset(ndata, 'Z', size);
    }
    MemCmd cmd = (mode==BaseTLB::Write) ? MemCmd::WriteReq :
        MemCmd::ReadReq;

    //virtual address request constructor (copied the data port request)
    //const int asid = 0;
    const Addr pc = tc->instAddr();
    RequestPtr req = std::make_shared<Request>(addr, size, 0,
        owner->VectorCacheMasterId, pc, tc->contextId());

    BaseCPU *cpu = tc->getCpuPtr();

    req->taskId(cpu->taskId());

    //start translation
    Tlb_Translation *translation = new Tlb_Translation(owner);

    tc->getDTBPtr()->translateTiming(req, tc, translation , mode);

    if (translation->fault == NoFault){
        PacketPtr pkt = new VectorPacket(req, cmd, req_id, channel);
        pkt->dataDynamic(ndata);

        if (!sendTimingReq(pkt)) {
            laCachePktQs[channel].push_back(pkt);
        }

        delete translation;
        return true;
    }else{
        translation->fault->invoke(tc, NULL);
        delete translation;
        return false;
        }
}

bool
VectorEngine::VectorMemPort::sendTimingReadReq(Addr addr, uint64_t size,
    ThreadContext *tc, uint64_t req_id, uint8_t channel)
{
    return startTranslation(addr, nullptr, size, BaseTLB::Read, tc, req_id,
        channel);
}

bool
VectorEngine::VectorMemPort::sendTimingWriteReq(Addr addr,
    uint8_t *data, uint64_t size, ThreadContext *tc, uint64_t req_id,
    uint8_t channel)
{
    return startTranslation(addr, data, size, BaseTLB::Write, tc, req_id,
        channel);
}

bool
VectorEngine::VectorMemPort::recvTimingResp(PacketPtr pkt)
{
    VectorPacketPtr la_pkt = dynamic_cast<VectorPacketPtr>(pkt);
    assert(la_pkt != nullptr);
    owner->recvTimingResp(la_pkt);
    return true;
}

void
VectorEngine::VectorMemPort::recvReqRetry()
{
    //TODO: must be a better way to figure out which channel specified the
    //      retry
    for (auto& laCachePktQ : laCachePktQs) {
        //assert(laCachePktQ.size());
        while (laCachePktQ.size() && sendTimingReq(laCachePktQ.front())) {
            // memory system takes ownership of packet
            laCachePktQ.pop_front();
        }
    }
}

VectorEngine::VectorRegPort::VectorRegPort(const std::string& name,
    VectorEngine* owner, uint64_t channel) :
    MasterPort(name, owner), owner(owner), channel(channel)
{
}

VectorEngine::VectorRegPort::~VectorRegPort()
{
}

bool
VectorEngine::VectorRegPort::sendTimingReadReq(Addr addr, uint64_t size,
    uint64_t req_id)
{
    //physical addressing
    RequestPtr req = std::make_shared<Request>
        (addr,size,0,owner->VectorRegMasterIds[channel]);
    PacketPtr pkt = new VectorPacket(req, MemCmd::ReadReq, req_id);

    //make data for cache to put data into
    uint8_t *ndata = new uint8_t[size];
    memset(ndata, 'Z', size);
    pkt->dataDynamic(ndata);

    if (!sendTimingReq(pkt)) {
        //delete pkt->req;
        delete ndata;
        delete pkt;
        return false;
    }
    return true;
}

// memcpy of data, so the caller can delete it when it pleases
bool
VectorEngine::VectorRegPort::sendTimingWriteReq(Addr addr, uint8_t *data,
    uint64_t size, uint64_t req_id)
{
    //physical addressing
    RequestPtr req = std::make_shared<Request>
        (addr,size,0,owner->VectorRegMasterIds[channel]);
    VectorPacketPtr pkt = new VectorPacket(req, MemCmd::WriteReq, req_id);

    //make copy of data here
    uint8_t *ndata = new uint8_t[size];
    memcpy(ndata, data, size);
    pkt->dataDynamic(ndata);

    if (!sendTimingReq(pkt)){
        //delete pkt->req;
        delete ndata;
        delete pkt;
        return false;
    }
    return true;
}

bool
VectorEngine::VectorRegPort::recvTimingResp(PacketPtr pkt)
{
    VectorPacketPtr vector_pkt = dynamic_cast<VectorPacketPtr>(pkt);
    assert(vector_pkt != nullptr);
    owner->recvTimingResp(vector_pkt);
    return true;
}

void
VectorEngine::VectorRegPort::recvReqRetry()
{
    //do nothing, caller will manually resend request
}

void
VectorEngine::recvTimingResp(VectorPacketPtr vector_pkt)
{
    //find the associated request in the queue and associate with vector_pkt
    assert(vector_PendingReqQ.size());
    bool found = false;
    for (Vector_ReqState * pending : vector_PendingReqQ) {
        if (pending->getReqId() == vector_pkt->reqId){
            pending->setPacket(vector_pkt);
            found = true;
            break;
        }
    }
    assert(found);

    //commit each request in order they were issued
    while (vector_PendingReqQ.size() &&
        vector_PendingReqQ.front()->isMatched())
    {
        Vector_ReqState * pending = vector_PendingReqQ.front();
        vector_PendingReqQ.pop_front();
        pending->executeCallback();
        //NOTE: pending deletes packet. packet calls delete[] on the data
        delete pending;
    }
}


bool
VectorEngine::writeVectorMem(Addr addr, uint8_t *data, uint32_t size,
    ThreadContext *tc, uint8_t channel, std::function<void(void)> callback)
{
    uint64_t id = (uniqueReqId++);
    Vector_ReqState *pending = new Vector_W_ReqState(id, callback);
    vector_PendingReqQ.push_back(pending);
    if (!vectormem_port.sendTimingWriteReq(addr, data, size, tc, id, channel)){
        delete vector_PendingReqQ.back();
        vector_PendingReqQ.pop_back();
        return false;
    }
    return true;
}

bool
VectorEngine::writeVectorReg(Addr addr, uint8_t *data,
    uint32_t size, uint8_t channel,
    std::function<void(void)> callback)
{
    //DPRINTF(VectorEngine, "writeVectorReg got %d bytes to write at %#x\n"
    //    ,size, addr);
    uint64_t id = (uniqueReqId++);
    Vector_ReqState *pending = new Vector_W_ReqState(id, callback);
    vector_PendingReqQ.push_back(pending);
    if (!VectorRegPorts[channel].sendTimingWriteReq(addr,data,size,id)) {
        delete vector_PendingReqQ.back();
        vector_PendingReqQ.pop_back();
        return false;
    }
    return true;
}

bool
VectorEngine::readVectorMem(Addr addr, uint32_t size,
    ThreadContext *tc, uint8_t channel,
    std::function<void(uint8_t*,uint8_t)> callback)
{
    uint64_t id = (uniqueReqId++);
    Vector_ReqState *pending = new Vector_R_ReqState(id, callback);
    vector_PendingReqQ.push_back(pending);
    if (!vectormem_port.sendTimingReadReq(addr, size, tc, id, channel)) {
        delete vector_PendingReqQ.back();
        vector_PendingReqQ.pop_back();
        return false;
    }
    return true;
}

bool
VectorEngine::readVectorReg(Addr addr, uint32_t size,
    uint8_t channel,
    std::function<void(uint8_t*,uint8_t)> callback)
{
    uint64_t id = (uniqueReqId++);
    Vector_ReqState *pending = new Vector_R_ReqState(id,callback);
    vector_PendingReqQ.push_back(pending);
    if (!VectorRegPorts[channel].sendTimingReadReq(addr,size,id)) {
        delete vector_PendingReqQ.back();
        vector_PendingReqQ.pop_back();
        return false;
    }
    return true;
}

VectorEngine *
VectorEngineParams::create()
{
    return new VectorEngine(this);
}