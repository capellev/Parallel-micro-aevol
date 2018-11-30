#pragma once


__device__
void dna_gpu_copy(char** dna, char** dna_next_gen, size_t dna_size, int indiv_id) {
    for (int i = 0; i < dna_size;i++)
        dna_next_gen[indiv_id][i] = dna[indiv_id][i];
}

__device__
void dna_gpu_set(char** dna, int indiv_id, int pos, char c) {
    dna[indiv_id][pos] = c;
}

__device__
void dna_gpu_do_switch(char** dna, int indiv_id, int pos) {
    if (dna[indiv_id][pos] == '0') dna[indiv_id][pos] = '1';
    else dna[indiv_id][pos] = '0';
}