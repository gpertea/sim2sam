#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>

#include "gff.h"
#include "GFaSeqGet.h"
#include "FastaTools.h"
#include "GFaSeqGet.h"
#include "arg_parse.h"

#define BAM_CMATCH  0
#define BAM_CINS    1
#define BAM_CDEL    2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD    6
#define BAM_CEQUAL  7
#define BAM_CDIFF   8
#define BAM_CIGAR_SHIFT 4
#define MAX_CIGARS  1024

class Position{
public:
    uint8_t num_elems = 0;
    std::string chr = "";
    char strand = '*';
    uint32_t start=0,locus=0;
    double abund = 0; // set abundance to a given position. This would require us to somehow tell which transcripts are shared by a given position
    int num_used = 0; // whn performing abundance allocation - how many times has the value been written - used for pid correction
    std::vector<std::string> transIDs; // for indexing, this includes all transcripts that contain this position. Otherwise can be any transcript which can describe the position to which it belongs; it does not participate in the equality computation
    std::vector<uint32_t> moves{}; // simplified CIGAR describing the intron-exon coverage of the given kmer
    bool revcmp=false;
    bool uniq_trans = true;
    int len = 0;

    Position(){
        num_elems = 0;
        chr="",strand='*';
        start=0,locus=0;
        len = 0;
        moves.resize(0); // simplified CIGAR describing the intron-exon coverage of the given kmer
        revcmp=false;
        abund=0;
        num_used=0;
        uniq_trans=true;
    };
    Position(std::string chr,char strand,uint32_t start,uint32_t locus,std::string transID,bool revcmp,int len){
        this->chr = chr;
        this->strand = strand;
        this->start = start;
        this->len = len;
        this->locus = locus;
        this->transIDs.push_back(transID);
        this->revcmp = revcmp;
    }
    ~Position() = default;

    void add_move(uint32_t move){this->moves.push_back(move);this->num_elems++;}
    void set_chr(std::string new_chr){this->chr=new_chr;this->num_elems++;}
    void set_strand(char new_strand){this->strand=new_strand;this->num_elems++;}
    void set_start(uint32_t new_start){this->start=new_start;this->num_elems++;}
    void set_locus(uint32_t locID){this->locus=locID;this->num_elems++;}
    void add_transID(std::string trans){
        this->transIDs.push_back(trans);
        this->num_elems++;
    }
    void set_rev(bool rev_tag){this->revcmp=rev_tag;this->num_elems++;}
    void set_uniq(bool uniq_tag){this->uniq_trans=uniq_tag;this->num_elems++;}

    static bool moves_eq(const std::vector<uint32_t>& m1,const std::vector<uint32_t>& m2) { // m1 must be smaller or equal to m2
        for(int i=0;i<m1.size();i++){
            if(m1[i] != m2[i]){
                return i == m1.size() - 1 && m1[i] <= m2[i];
            }
        }
        return true;
    }

    std::string get_strg() const {
        // remove duplicates from the vector of transIDs
        std::set<std::string> tmp_transIDs;
        for(auto t : transIDs){
            tmp_transIDs.insert(t);
        }

        std::string res;
        res.append(std::to_string(this->revcmp)); // the first character always indicates whether reverse complemented or not
        res.append(std::to_string(this->uniq_trans)); // the second character always indicates on whether the multimapper should be treated as having come from a sufficiently percent unique locus or not
        for(auto t : tmp_transIDs){
            res.append(t);
            res += '^';
        }
        res.pop_back(); // remove the last ^
        res += ">";
        res.append(std::to_string(this->locus));
        res += "@";
        res.append(this->chr);
        res += this->strand;
        res.append(std::to_string(this->start));
        res += ':';
        for(auto &mit : this->moves){
            res.append(std::to_string(mit));
            res += ' ';
        }
        res.pop_back(); // removes the last whitespace
        return res;
    }

    bool operator==(const Position& m) const{
        return this->chr==m.chr &&
               this->strand==m.strand &&
               this->start==m.start &&
               this->moves==m.moves;
    }

    bool operator>(const Position& m) const{
        return this->locus > m.locus;
    }

    bool operator<(const Position& m) const{
        if(this->chr<m.chr){
            return true;
        }
        if(this->chr>m.chr){
            return false;
        }
        if(this->strand<m.strand){
            return true;
        }
        if(this->strand>m.strand){
            return false;
        }
        if(this->start<m.start){
            return true;
        }
        if(this->start>m.start){
            return false;
        }
        if(this->moves<m.moves){
            return true;
        }
    }

    bool lt_noStrand(const Position& rhs) const{
        if(this->chr<rhs.chr){
            return true;
        }
        if(this->chr>rhs.chr){
            return false;
        }
        if(this->start<rhs.start){
            return true;
        }
        if(this->start>rhs.start){
            return false;
        }
        if(this->moves<rhs.moves){
            return true;
        }
    }

    void clear(){
        this->moves.clear();
        this->transIDs.clear();
        this->num_elems = 0;
        this->revcmp = false;
        this->abund=0;
        this->num_used=0;
        this->uniq_trans=true;
    }

    uint8_t size(){return this->num_elems;}

    Position(const Position &p2) {
        num_elems = p2.num_elems;
        chr = p2.chr;
        strand = p2.strand;
        start = p2.start;
        locus = p2.locus;
        abund = p2.abund;
        num_used=p2.num_used;
        for(auto t : p2.transIDs){
            this->add_transID(t);
        }
        for(auto m : p2.moves){
            this->add_move(m);
        }
        this->revcmp = p2.revcmp;
        this->uniq_trans = p2.uniq_trans;
    }
};

// this function generates alignment header
void get_header(std::string &header,std::string index_fname,std::string cl){
    header += "@HD\tVN:1.0\tSO:unsorted\n";
    // read the index and extract information
    std::ifstream index_stream;
    std::stringstream *linestream;
    index_stream.open(index_fname.c_str(),std::ios::in);
    if (!index_stream.good()){
        std::cerr<<"@ERROR::Couldn't open FASTA index: "<<index_fname<<std::endl;
        exit(1);
    }
    std::ios::sync_with_stdio(false);

    std::string chrID,chrLen;
    std::string aline,col;
    while (std::getline(index_stream,aline)) {
        // given a line we need to extract the name and the TPM
        linestream = new std::stringstream(aline);
        std::getline(*linestream,col,'\t');
        chrID = col;
        // now need to get the abundance
        std::getline(*linestream,col,'\t'); // skip second column
        chrLen = col;
        header += "@SQ\tSN:" + chrID + "\tLN:" + chrLen + "\n";
        delete linestream;
    }
    index_stream.close();
    // PG is added using the arguments from the execution of trans2genome during the conversion
    header += "@PG\tID:sim2sam\tPN:sim2sam:0.0\tCL:\"" + cl + "\"";
}

class Read{
public:
    Read():name(""),tid(""),start(0),end(0),strand(-1){}
    ~Read() = default;
    void parse_read_poly(std::string &read){
        // get read number which is first
        std::size_t prev_delim = 0;
        std::size_t delim = read.find("/");
        this->name = read.substr(0, delim);

        // get transcript id
        prev_delim = delim;
        delim = read.find(";", prev_delim + 1);
        this->tid = read.substr(prev_delim + 1, (delim - prev_delim) - 1);

        // get start position on the transcript
        prev_delim = delim;
        delim = read.find(":", prev_delim + 1);
        prev_delim = delim;
        delim = read.find("-", prev_delim + 1);
        this->start = std::atoi(read.substr(prev_delim + 1, (delim - prev_delim) - 1).c_str()) - 1;

        // get end position on the transcript
        prev_delim = delim;
        delim = read.find(";", prev_delim + 1);
        this->end = std::atoi(read.substr(prev_delim + 1, (delim - prev_delim) - 1).c_str());
    }
    void parse_read_rsem(std::string &read,std::vector<std::string> &rsemi){
        std::size_t prev_delim = 0;
        std::size_t delim = read.find('_');
        this->name = read.substr(0,delim);

        prev_delim = delim;
        delim = read.find('_',prev_delim+1);
        this->strand = std::atoi(read.substr(prev_delim+1,(delim-prev_delim)-1).c_str());

        prev_delim = delim;
        delim = read.find('_',prev_delim+1);
        int tmp_id = std::atoi(read.substr(prev_delim+1,(delim-prev_delim)-1).c_str());
        this->tid = rsemi[tmp_id-1];

        prev_delim = delim;
        prev_delim = delim;
        delim = read.find('_',prev_delim+1);
        this->start = std::atoi(read.substr(prev_delim+1,(delim-prev_delim)-1).c_str());
    }
    void clear(){
        name = "";
        tid = "";
        start = 0;
        end = 0;
        strand = -1;
    }

    std::string name,tid;
    int start,end,type,strand;
};

bool get_read_start(GList<GffExon>& exon_list,int tstart,int &gstart,int &exon_idx){
    if(tstart==-1){ // conforming to SAM specification where 0 pos PNEXT means it is not set (for bowtie means reads are not paired)
        gstart = -1;
        return true;
    }
    const GSeg* cur_exon;
    size_t cur_intron_dist=0;
    size_t trans_start=exon_list[0]->start;
    int trans_offset=0;
    for(int i=0;i<exon_list.Count();++i){
        cur_exon= exon_list[i];
        trans_offset=trans_start + cur_intron_dist;

        if (tstart>=cur_exon->start - trans_offset && tstart <=cur_exon->end - trans_offset){
            gstart = tstart + trans_start + cur_intron_dist;
            exon_idx=i;
            return true;
        }
        else{
            if (i+1<exon_list.Count()){
                cur_intron_dist += exon_list[i+1]->start - cur_exon->end-1;
            }
            else{
                return false;
            }
        }
    }
    return false;
}

bool get_read_start_neg(GList<GffExon>& exon_list,int tstart,int &gstart,int &exon_idx){
    if(tstart==-1){ // conforming to SAM specification where 0 pos PNEXT means it is not set (for bowtie means reads are not paired)
        gstart = -1;
        return true;
    }
    const GSeg* cur_exon;
    size_t cur_intron_dist=0;
    size_t trans_end=exon_list[exon_idx]->end;
    int trans_offset=0;
    for(int i=exon_idx;i>=0;--i){
        cur_exon = exon_list[i];
        trans_offset = trans_end - cur_intron_dist;

        if (tstart>=-(cur_exon->end - trans_offset) && tstart <=cur_exon->end - trans_offset){
            gstart = trans_end - cur_intron_dist - tstart;
            exon_idx=i;
            return true;
        }
        else{
            if (i-1>=0){
                cur_intron_dist += (exon_list[i-1]->end-1) - cur_exon->start;
            }
            else{
                return false;
            }
        }
    }
    return false;
}

int get_cigar(int exon_i,GSeg *next_exon,GList<GffExon>& exon_list,std::string &cigars,int read_start,int readlen){
    int cur_pos = read_start;
    int remaining_length = readlen;
    int miss_length = 0,match_length = 0;

    for (;exon_i<exon_list.Count();++exon_i){
        GffExon* cur_exon=exon_list[exon_i];
        if (cur_pos>=(int)cur_exon->start && cur_pos+remaining_length-1<=(int)cur_exon->end){
//            cigars[num_cigars] = BAM_CMATCH | (remaining_length <<BAM_CIGAR_SHIFT);
//            ++num_cigars;
            cigars += std::to_string(remaining_length) + "M";
            cur_pos+=remaining_length;
            break;
        }
        else if (cur_pos >= (int)cur_exon->start && cur_pos+remaining_length-1>(int)cur_exon->end){
            match_length=(int)cur_exon->end-cur_pos+1;
            if (match_length>0){
//                cigars[num_cigars]=BAM_CMATCH | (match_length << BAM_CIGAR_SHIFT);
//                ++num_cigars;
                cigars += std::to_string(match_length) + "M";
            }
            if (exon_i+1>=exon_list.Count()){
                return 0;
            }
            else{
                next_exon=exon_list[exon_i+1];
            }
            miss_length=next_exon->start-cur_exon->end-1;

//            cigars[num_cigars]=BAM_CREF_SKIP | (miss_length <<BAM_CIGAR_SHIFT);
//            ++num_cigars;
            cigars += std::to_string(miss_length) + "N";

            cur_pos+=match_length+miss_length;

            remaining_length-=match_length;
            assert(cur_pos == (int)next_exon->start);
        }
    }
}

void process_read(Position &cur_pos,std::string &cigars,GffObj *p_gffObj,int tstart,int read_len,char actual_strand){
    GList<GffExon>& exon_list = p_gffObj->exons; // get exons

    GSeg *next_exon=nullptr;
    int exon_i = (actual_strand == '+') ? 0 : exon_list.Count()-1;
    int32_t gstart=0;

    // first find the genomic read start
    if(actual_strand == '-'){
        // compute total transcript length
        int tlen = 0;
        for(int i=0;i<p_gffObj->exons.Count();i++){
            tlen+=p_gffObj->exons[i]->len();
        }
        tstart = tlen - (tstart+read_len);
    }
    bool ret_val = get_read_start(exon_list,tstart,gstart,exon_i);
    if(!ret_val){
        std::cerr<<"@ERROR::Can not get genomic read start"<<std::endl;
        exit(1);
    }

    cur_pos.set_chr(std::string(p_gffObj->getRefName()));
    cur_pos.set_start(gstart);
    cur_pos.set_strand(actual_strand);
    cur_pos.add_transID(std::string(p_gffObj->getID()));

    // secondly build a new cigar string
    ret_val = get_cigar(exon_i,next_exon,exon_list,cigars,gstart,read_len);
    if (!ret_val) {
        std::cerr << "@ERROR::Can not create a new cigar string for the single read from process_read" << std::endl;
        exit(1);
    }
}

enum Opt {GFF = 'g',
    SINGLE    = 's',
    OUTPUT    = 'o',
    INDEX     = 'i',
    RSEM_MAP  = 'r',
    TYPE      = 't'
};


int main(int argc, char** argv) {
    ArgParse args("sim2sam");
    args.add_string(Opt::GFF,"gff","","annotation from which reads are simulated",true);
    args.add_string(Opt::SINGLE,"single","","single-end reads",true);
    args.add_string(Opt::OUTPUT,"output","","output SAM alignment",true);
    args.add_string(Opt::INDEX,"index","","Fasta index of the genome from which reads where simulated",true);
    args.add_string(Opt::RSEM_MAP,"rsemi","","Base name of the RSEM-prepared reference. Required if the RSEM mode is enabled",false);
    args.add_string(Opt::TYPE,"type","","Type of the simulator used. Options are rsem or polyester",true);

    if(argc <= 1 || strcmp(argv[1],"--help")==0){
        std::cerr<<args.get_help()<<std::endl;
        exit(1);
    }

    if(args.get_string(Opt::TYPE)=="rsem"){
        if(!args.is_set(Opt::RSEM_MAP)){
            std::cerr<<"To parse RSEM simulation --rsemi must be specified"<<std::endl;
            exit(1);
        }
    }

    args.parse_args(argc,argv);

    // first create the execution string
    std::string cl="sim2sam ";
    for (int i=0;i<argc;i++){
        if(i==0){
            cl+=argv[i];
        }
        else{
            cl+=" ";
            cl+=argv[i];
        }
    }

    std::ofstream out_al_fp(args.get_string(Opt::OUTPUT));

    std::string header = "";
    get_header(header,args.get_string(Opt::INDEX),cl);
    // write header
    out_al_fp << header << std::endl;

    // load the GFF
    FILE* gff_file = fopen(args.get_string(Opt::GFF).c_str(), "r");
    if (gff_file == nullptr)
    {
        std::cerr << "@ERROR::Couldn't open annotation: " << args.get_string(Opt::GFF)<< std::endl;
        exit(1);
    }
    GffReader gffReader;
    gffReader.init(gff_file,true);
    gffReader.readAll();
    GffObj *p_gffObj;
    // build a map of all transcript IDs to indices within gfflst
    std::unordered_map<std::string,int> idm;
    std::pair<std::unordered_map<std::string,int>::iterator,bool> idm_it;
    for (int i = 0; i < gffReader.gflst.Count(); ++i){
        p_gffObj = gffReader.gflst.Get(i);
        if (p_gffObj->isDiscarded() || p_gffObj->exons.Count()==0) continue;
        idm_it = idm.insert(std::make_pair(p_gffObj->getID(),i));
        if(!idm_it.second){
            std::cerr<<"duplicate transcript IDs detected: "<<p_gffObj->getID()<<" at index "<<i<<std::endl;
            exit(-1);
        }
    }
    std::cerr<<"number of transcripts in the index is: "<<idm.size()<<std::endl;

    // if type is RSEM then need to preload the index
    std::vector<std::string> rsemi;
    if(args.get_string(Opt::TYPE)=="rsem"){
        std::vector<int> tid_tmps;

        std::ifstream rsem_idx_fp(args.get_string(Opt::RSEM_MAP));
        std::string line;
        std::getline(rsem_idx_fp,line); // skip header
        int line_no = 0;
        while (std::getline(rsem_idx_fp, line)){
            std::istringstream iss(line);
            std::string tid;
            std::getline( iss >> std::skipws,tid,'\t');
            line_no++;
            rsemi.push_back(tid);
        }
        rsem_idx_fp.close();
    }

    enum TP {RSEM = 0,
             POLY = 1
    };
    int type;
    if(args.get_string(Opt::TYPE)=="rsem"){
        type = TP::RSEM;
    }
    else if(args.get_string(Opt::TYPE)=="polyester"){
        type = TP::POLY;
    }
    else{
        std::cerr<<"unknown type specified"<<std::endl;
        exit(1);
    }

    // now need to iterate over the reads and write out SAM records
    // TODO: need same implementation for paired reads
    FastaReader fastaReader(args.get_string(Opt::SINGLE));
    FastaRecord read_str;
    Read read;
    int read_len;
    char strand;
    int count = 0;
    while (fastaReader.good()) {
        fastaReader.next(read_str);
        if(type==TP::RSEM){
            read.parse_read_rsem(read_str.id_,rsemi);
        }
        else if(type==TP::POLY){
            read.parse_read_poly(read_str.id_);
        }
        else{
            std::cerr<<"unsupported type"<<std::endl;
            exit(-1);
        }

        // get transcript from gffReader
        idm_it.first = idm.find(read.tid); // get index within the gfflst
        if(idm_it.first == idm.end()){
            std::cerr<<"invalid transcript id found in the reads: "<<read.tid<<std::endl;
            exit(-1);
        }
        p_gffObj = gffReader.gflst.Get(idm_it.first->second); // get the actual transcript object
        strand = p_gffObj->strand;
        if(type==TP::RSEM){ // strand is also set within the read and neads to be taken into account
            if(strand=='+' && read.strand){
                strand = '-';
            }
            else if(strand=='+' && !read.strand){
                strand = '+';
            }
            else if(strand=='-' && read.strand){
                strand = '+';
            }
            else if(strand=='-' && !read.strand){
                strand = '-';
            }
            else{
                std::cerr<<"unknown strand"<<std::endl;
                exit(-1);
            }
        }
        // not to get the genomic position and the cigar string
        Position cur_pos;
        std::string cigars = "";
        read_len = read_str.seq_.size();
        process_read(cur_pos,cigars,p_gffObj,read.start,read_len,strand);

        // get flag and reverse sequence
        int flag = 0;
        char seq[read_str.seq_.size() + 1];
        strcpy(seq, read_str.seq_.c_str());
        if(cur_pos.strand == '-'){
            flag = 16;
            reverseComplement(seq,read_str.seq_.size());
        }

        // now to write the record to the file
        out_al_fp << read.name << "\t"
                  << flag << "\t"
                  << std::string(p_gffObj->getRefName()) << "\t"
                  << cur_pos.start << "\t"
                  << 60 << "\t"
                  << cigars <<"\t"
                  << "*" << "\t"
                  << "0" << "\t"
                  << "0" <<"\t"
                  << std::string(seq) << "\t"
                  << std::string(read_len,'I') << "\t"
                  << "AS:i:0" << "\t"
                  << "XN:i:0" << "\t"
                  << "XM:i:0" << "\t"
                  << "XO:i:0" << "\t"
                  << "XG:i:0" << "\t"
                  << "NM:i:0" << "\t"
                  << "MD:Z:100" << "\t"
                  << "YT:Z:UU" << "\t"
                  << "NH:i:1" << "\t"
                  << "XS:A:" << cur_pos.strand
                  << std::endl;

        read.clear();
    }

    out_al_fp.close();
    return 0;
}