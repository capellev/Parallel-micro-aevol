//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"
#include "ExpManager.h"

Dna::Dna(const Dna& clone) : seq_(clone.seq_) {
}

Dna::Dna(int length, Threefry::Gen& rng) : seq_(length) {
  // Generate a random genome
  for (int32_t i = 0; i < length; i++) {
    seq_[i] = '0' + rng.random(NB_BASE);
  }
}

Dna::Dna(char* genome, int length) : seq_(length) {
  strcpy(seq_.data(), genome);
}

Dna::Dna(int length) : seq_(length) {
}

int Dna::length() const {
  return seq_.size();
}

void Dna::save(gzFile backup_file) {
  int dna_length = length();
  gzwrite(backup_file, &dna_length, sizeof(dna_length));
  gzwrite(backup_file, seq_.data(), seq_.size() * sizeof(seq_[0]));
}

void Dna::load(gzFile backup_file) {
  int dna_length;
  gzread(backup_file,&dna_length,sizeof(dna_length));

  gzread(backup_file, seq_.data(), sizeof(*seq_.data()));
}

void Dna::set(int pos, char c) {
  seq_[pos] = c;
}

void Dna::do_switch(int pos) {
  if (seq_[pos] == '0') seq_[pos] = '1';
  else seq_[pos] = '0';
}


int Dna::promoter_at(int pos) {
  int prom_dist[22];

  for (int motif_id = 0; motif_id < 22; motif_id++) {
    // Searching for the promoter
    prom_dist[motif_id] =
        PROM_SEQ[motif_id] ==
        seq_[
            pos + motif_id >= seq_.size() ? pos +
                                            motif_id -
                                            seq_.size()
                                          : pos +
                                            motif_id]
        ? 0
        : 1;

  }


  // Computing if a promoter exists at that position
  int dist_lead = prom_dist[0] +
                  prom_dist[1] +
                  prom_dist[2] +
                  prom_dist[3] +
                  prom_dist[4] +
                  prom_dist[5] +
                  prom_dist[6] +
                  prom_dist[7] +
                  prom_dist[8] +
                  prom_dist[9] +
                  prom_dist[10] +
                  prom_dist[11] +
                  prom_dist[12] +
                  prom_dist[13] +
                  prom_dist[14] +
                  prom_dist[15] +
                  prom_dist[16] +
                  prom_dist[17] +
                  prom_dist[18] +
                  prom_dist[19] +
                  prom_dist[20] +
                  prom_dist[21];

  return dist_lead;
}

int Dna::terminator_at(int pos) {
  int term_dist[4];
  for (int motif_id = 0; motif_id < 4; motif_id++) {

    // Search for the terminators
    term_dist[motif_id] =
        seq_[
            pos + motif_id >= seq_.size() ? pos +
                                            motif_id -
                                            seq_.size() :
            pos + motif_id] !=
        seq_[
            pos - motif_id + 10 >= seq_.size() ?
            pos - motif_id + 10 - seq_.size() :
            pos -
            motif_id +
            10] ? 1
                : 0;
  }
  int dist_term_lead = term_dist[0] +
                       term_dist[1] +
                       term_dist[2] +
                       term_dist[3];

  return dist_term_lead;
}

bool Dna::shine_dal_start(int pos) {
  bool start = false;
  int t_pos, k_t;

  for (int k = 0; k < 9; k++) {
    k_t = k >= 6 ? k + 4 : k;
    t_pos = pos + k_t >= seq_.size() ? pos + k_t -
                                       seq_.size()
                                     : pos + k_t;

    if (seq_[t_pos] ==
        SHINE_DAL_SEQ[k]) {
      start = true;
    } else {
      start = false;
      break;
    }
  }

  return start;
}

bool Dna::protein_stop(int pos) {
  bool is_protein;
  int t_k;

  for (int k = 0; k < 3; k++) {
    t_k = pos + k >= seq_.size() ?
          pos - seq_.size() + k :
          pos + k;

    if (seq_[t_k] ==
        PROTEIN_END[k]) {
      is_protein = true;
    } else {
      is_protein = false;
      break;
    }
  }

  return is_protein;
}

int Dna::codon_at(int pos) {
  int value = 0;

  int t_pos;

  for (int i = 0; i < 3; i++) {
    t_pos =
        pos + i >= seq_.size() ? pos + i -
                                 seq_.size()
                               : pos + i;
    if (seq_[t_pos] ==
        '1')
      value += 1 << (CODON_SIZE - i - 1);
  }

  return value;
}