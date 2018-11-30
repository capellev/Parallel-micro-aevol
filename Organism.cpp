// ***************************************************************************************************************
//
//          Mini-Aevol is a reduced version of Aevol -- An in silico experimental evolution platform
//
// ***************************************************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <https://gitlab.inria.fr/rouzaudc/mini-aevol>
// Web: https://gitlab.inria.fr/rouzaudc/mini-aevol
// E-mail: See <jonathan.rouzaud-cornabas@inria.fr>
// Original Authors : Jonathan Rouzaud-Cornabas
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ***************************************************************************************************************



#include <cstring>
#include "Organism.h"
#include "ExpManager.h"

#include <iostream>
using namespace std;

/**
 * Constructor to generate a random organism (i.e. an organism with a random DNA)
 *
 * @param exp_m : Related ExpManager object
 * @param length : Length of the generated random DNA
 * @param indiv_id : Unique Identification Number
 */
Organism::Organism(ExpManager* exp_m, int length, int indiv_id) {
    exp_m_ = exp_m;

    count_prom = 0;
    rna_count_ = 0;

    auto rng = exp_m->rng_->gen(indiv_id, Threefry::MUTATION);
    dna_ = new Dna(length, rng);
    parent_length_ = length;
    indiv_id_ = indiv_id;
}

/**
 * Create an organism with a given genome
 *
 * @param exp_m : Related ExpManager object
 * @param genome : Genome to assign to the organism
 * @param indiv_id : Unique Identification Number
 */
Organism::Organism(ExpManager *exp_m, char* genome, int indiv_id) {
    exp_m_ = exp_m;

    count_prom = 0;
    rna_count_ = 0;

    dna_ = new Dna(genome,strlen(genome));
    parent_length_ = strlen(genome);
    indiv_id_ = indiv_id;

}

/**
 * Constructor to create a clone of a given Organism
 *
 * @param exp_m : Related ExpManager object
 * @param clone : The organism to clone
 */
Organism::Organism(ExpManager *exp_m, std::shared_ptr<Organism> clone) {
    exp_m_ = exp_m;

    count_prom = 0;
    rna_count_ = 0;

  parent_length_ = clone->parent_length_;
    dna_ = new Dna(*(clone->dna_));

    for (const auto& prom : clone->promoters) {
        if (prom.second != nullptr) {
            auto prom_copy = new Promoter(prom.second->pos, prom.second->error);
            promoters[count_prom] = prom_copy;

            prom_pos[prom_copy->pos] = count_prom;

            count_prom++;
        }
    }
}

/**
 * Create an Organism from a backup/checkpointing file
 *
 * @param exp_m : Related ExpManager object
 * @param backup_file : gzFile to read from
 */
Organism::Organism(ExpManager *exp_m, gzFile backup_file) {
    exp_m_ = exp_m;

    count_prom = 0;
    rna_count_ = 0;

    load(backup_file);

}

/**
 * Destructor of an organism
 */
Organism::~Organism() {
    for (auto prom : promoters) {
        delete(prom.second);
    }
    promoters.clear();
    prom_pos.clear();
    
    for (auto rna : rnas) {
        delete(rna);
    }
    rnas.clear();
    
    for (auto prot : proteins) {
        delete(prot);
    }
    proteins.clear();
    
    terminators.clear();
    
    delete dna_;
    
}

/**
 * Save the organism to backup/checkpointing file
 *
 * @param backup_file : where to the save the organism
 */
void Organism::save(gzFile backup_file) {
    gzwrite(backup_file,&indiv_id_,sizeof(indiv_id_));
    gzwrite(backup_file,&parent_id_,sizeof(parent_id_));
    gzwrite(backup_file,&global_id,sizeof(global_id));
;
    gzwrite(backup_file, &parent_length_, sizeof(parent_length_));

  dna_->save(backup_file);
}

/**
 * Load the organism from backup/checkpointing file
 *
 * @param backup_file : from where restore the organism
 */
void Organism::load(gzFile backup_file) {
    gzread(backup_file,&indiv_id_,sizeof(indiv_id_));
    gzread(backup_file,&parent_id_,sizeof(parent_id_));
    gzread(backup_file,&global_id,sizeof(global_id));

    int parent_length;
    gzread(backup_file,&parent_length,sizeof(parent_length));

  dna_ = new Dna();
  dna_->load(backup_file);
}

/**
 * Reset the stats variable (used at the beginning of a generation when an organism is a perfect clone of its parent)
 */
void Organism::reset_stats() {
    nb_genes_activ = 0;
    nb_genes_inhib = 0;
    nb_func_genes = 0;
    nb_non_func_genes = 0;
    nb_coding_RNAs = 0;
    nb_non_coding_RNAs = 0;
}

/**
 * Replace the sequence of the DNA of the organism at a given position by a given sequence
 *
 * @param pos : where to replace the DNA by the given sequence
 * @param seq : the sequence itself
 * @param seq_length : length of the sequence
 */
void Organism::replace(int pos, char* seq, int seq_length) {
// Invert the sequence between positions 'first' and 'last'
    // Check pos value
    assert(pos >= 0 && pos < dna_->length());

    // If the sequence's length was not provided, compute it
    if (seq_length == -1) {
        seq_length = strlen(seq);
    }

    // Check that the sequence is contiguous
    assert(pos + seq_length <= dna_->length());

    // Perform the replacement
    memcpy(&dna_[pos], seq, seq_length * sizeof(char));
}

/**
 * Switch the DNA base-pair at a given position
 *
 * @param pos : the position where to switch the base-pair
 * @return
 */
bool Organism::do_switch(int pos) {
    dna_->do_switch(pos);

    // Remove promoters containing the switched base
    remove_promoters_around(pos, mod(pos + 1, length()));

    // Look for potential new promoters containing the switched base
    if (length() >= PROM_SIZE)
        look_for_new_promoters_around(pos, mod(pos + 1, length()));

    return true;
}

/**
 * Apply all the mutation events of the organism on its DNA
 */
void Organism::apply_mutations() {
    MutationEvent *repl;

    do {
        repl = exp_m_->
                dna_mutator_array_[indiv_id_]->generate_next_mutation(length());

        if (repl != nullptr) {
            switch (repl->type()) {
                case DO_SWITCH:
                    do_switch(repl->pos_1());
                    nb_swi_++;
                    nb_mut_++;
                    break;
            }
        }

    } while (exp_m_->dna_mutator_array_[indiv_id_]->mutation_available() > 0);

}

/**
Optimize promoters search
 **/


void Organism::remove_promoters_around(int32_t pos) {
    if (dna_->length() >= PROM_SIZE) {
        remove_promoters_starting_between(mod(pos - PROM_SIZE + 1,
                                                             dna_->length()),
                                                  pos);
    }
    else {
        remove_all_promoters();
    }
}

void Organism::remove_promoters_around(int32_t pos_1, int32_t pos_2) {
    if (mod(pos_1 - pos_2, dna_->length()) >= PROM_SIZE) {
        remove_promoters_starting_between(mod(pos_1 - PROM_SIZE + 1,
                                                             dna_->length()),
                                                  pos_2);
    }
    else {
        remove_all_promoters();
    }
}


void Organism::move_all_promoters_after(int32_t pos, int32_t delta_pos) {
    std::map<int32_t,int32_t> tmp_prom;

    for (auto it = prom_pos.lower_bound(pos), nextit=it;
         it != prom_pos.end();
         it = nextit) {

        int32_t new_pos = mod(it->first + delta_pos, dna_->length());
        int32_t prom_idx = it->second;

        promoters[it->second]->pos = new_pos;
        nextit = next(it);

        if (tmp_prom.find(new_pos) == tmp_prom.end()) {
            tmp_prom[new_pos] = prom_idx;
        } else {
            promoters.erase(it->second);
        }
        prom_pos.erase(it);

    }

    for (auto to_insert : tmp_prom) {
        if (prom_pos.find(to_insert.first) == prom_pos.end()) {
            prom_pos[to_insert.first] = to_insert.second;
        } else {
            promoters.erase(to_insert.second);
        }
    }
}

void Organism::look_for_new_promoters_around(int32_t pos_1, int32_t pos_2) {
    if (dna_->length() >= PROM_SIZE) {
        look_for_new_promoters_starting_between(
                mod(pos_1 - PROM_SIZE + 1,
                           dna_->length()), pos_2);
    }
}

void Organism::look_for_new_promoters_around(int32_t pos) {
    if (dna_->length() >= PROM_SIZE) {
        look_for_new_promoters_starting_between(
                mod(pos - PROM_SIZE + 1, dna_->length()),
                pos);
    }
}

void Organism::insert_promoters_at(std::list<Promoter*>& promoters_to_insert, int32_t pos) {

    if (promoters_to_insert.size() <= 0) {
        return;
    }

    // Insert the promoters in the individual's RNA list
    for (auto &to_insert: promoters_to_insert) {
        int prev_pos = to_insert->pos;
        // Update promoter position
        to_insert->pos = mod(to_insert->pos + pos, dna_->length());
        if (prom_pos.find(to_insert->pos) == prom_pos.end()) {
            int prom_idx = count_prom;
            count_prom = count_prom + 1;

            promoters[prom_idx] = to_insert;
            prom_pos[to_insert->pos] = prom_idx;
        }


    }
}


void Organism::duplicate_promoters_included_in(int32_t pos_1,
                                                           int32_t pos_2,
                                                           std::list<Promoter*>& duplicated_promoters) {
    // 1) Get promoters to be duplicated
    std::list<Promoter *> retrieved_promoters = {};

    promoters_included_in(pos_1, pos_2, retrieved_promoters);

    // 2) Set RNAs' position as their position on the duplicated segment
    for (auto &prom : retrieved_promoters) {
        // Make a copy of current RNA inside container
        duplicated_promoters.push_back(new Promoter(prom));

        // Set RNA's position as it's position on the duplicated segment
        duplicated_promoters.back()->pos = mod(duplicated_promoters.back()->pos - pos_1,
                                               dna_->length());
    }
}

void Organism::extract_promoters_included_in(int32_t pos_1,
                                                         int32_t pos_2,
                                                         std::list<Promoter*>& extracted_promoters) {
    if (pos_2 - pos_1 < PROM_SIZE) {
        return;
    }

    extract_promoters_starting_between(pos_1, pos_2 - PROM_SIZE + 1,
                                               extracted_promoters);
}

void Organism::insert_promoters(std::list<Promoter*>& promoters_to_insert) {
        if (promoters_to_insert.size() <= 0) {
            return;
        }

        // Insert the promoters in the individual's RNA list
        for (auto& to_insert: promoters_to_insert) {
                if (prom_pos.find(to_insert->pos) == prom_pos.end()) {

                    int prom_idx = count_prom;
                    count_prom = count_prom + 1;

                    promoters[prom_idx] = to_insert;
                    prom_pos[to_insert->pos] = prom_idx;
                }

        }

}

void Organism::remove_all_promoters() {
    prom_pos.clear();

    for (auto it = promoters.begin(),
                 nextit = it;
         it != promoters.end();
         it = nextit) {
        delete it->second;
        nextit = next(it);
        promoters.erase(it);
    }

    promoters.clear();
    count_prom = 0;
}

/** LEADING promoters **/
/** REMOVE **/
void Organism::remove_promoters_starting_between(int32_t pos_1, int32_t pos_2) {
    if (pos_1 > pos_2) {
        remove_promoters_starting_after(pos_1);
        remove_promoters_starting_before(pos_2);
    }
    else {
        // STL Warning: don't erase the current iterator in the for-loop!
        for (auto it = prom_pos.lower_bound(pos_1),
                     nextit = it;
             it != prom_pos.end() and it->first < pos_2;
             it = nextit) {

            int pidx = it->second;
            auto it_p = promoters[pidx];

            delete it_p;
            promoters.erase(pidx);
            nextit = next(it);
            prom_pos.erase(it);
        }
    }
}

void Organism::remove_promoters_starting_after(int32_t pos) {
    auto init_it = prom_pos.lower_bound(pos);
    if (init_it == prom_pos.end())
        return;

    for (auto it = init_it,
                 nextit = it;
         it != prom_pos.end();
         it = nextit) {
        delete promoters[it->second];
        promoters.erase(it->second);
        nextit = next(it);
        prom_pos.erase(it);
    }
}

void Organism::remove_promoters_starting_before(int32_t pos) {
// Delete RNAs until we reach pos (or we reach the end of the list)
    for (auto it = prom_pos.begin(),
                 nextit = it;
         it != prom_pos.end() and it->first < pos;
         it = nextit) {
        delete promoters[it->second];
        promoters.erase(it->second);
        nextit = next(it);
        prom_pos.erase(it);
    }
}


/** LOOK **/
void Organism::locate_promoters() {
    look_for_new_promoters_starting_between(0,dna_->length());
}

void Organism::look_for_new_promoters_starting_between(int32_t pos_1,int32_t pos_2) {
    // When pos_1 > pos_2, we will perform the search in 2 steps.
    // As positions  0 and dna_->length() are equivalent, it's preferable to
    // keep 0 for pos_1 and dna_->length() for pos_2.

    if (pos_1 >= pos_2) {
        look_for_new_promoters_starting_after(pos_1);
        look_for_new_promoters_starting_before(pos_2);
        return;
    }
    // Hamming distance of the sequence from the promoter consensus

    for (int32_t i = pos_1; i < pos_2; i++) {
        int8_t dist = dna_->promoter_at(i);

        if (dist <= 4) {
            if (prom_pos.find(i) == prom_pos.end()) {
                Promoter* nprom = new Promoter(i, dist);
                int prom_idx = count_prom;
                count_prom = count_prom + 1;

                promoters[prom_idx] = nprom;
                prom_pos[i] = prom_idx;
            }
        }
    }
}

void Organism::look_for_new_promoters_starting_after(int32_t pos) {
    for (int32_t i = pos; i < dna_->length(); i++) {
        int dist = dna_->promoter_at(i);

        if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
            if (prom_pos.find(i) == prom_pos.end()) {
                Promoter* nprom = new Promoter(i, dist);
                int prom_idx = count_prom;
                count_prom = count_prom + 1;

                promoters[prom_idx] = nprom;
                prom_pos[i] = prom_idx;
            }
        }
    }
}

void Organism::look_for_new_promoters_starting_before(int32_t pos) {
    // Hamming distance of the sequence from the promoter consensus

    for (int32_t i = 0; i < pos; i++) {

        int dist = dna_->promoter_at(i);

        if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
            if (prom_pos.find(i) == prom_pos.end()) {
                Promoter* nprom = new Promoter(i, dist);
                int prom_idx = count_prom;
                count_prom = count_prom + 1;

                promoters[prom_idx] = nprom;
                prom_pos[i] = prom_idx;
                }
            }
        }

}

/** EXTRACT **/
void Organism::extract_promoters_starting_between(int32_t pos_1,
                                                                      int32_t pos_2, std::list<Promoter*>& extracted_promoters) {
    if (pos_2 < pos_1) {

        auto first = prom_pos.lower_bound(pos_1);

        if (first == prom_pos.end() or first->first >= pos_2) {
            return;
        }

        // Extract the promoters (remove them from the individual's list and put them in extracted_promoters)

        for (auto it = first;
             it != prom_pos.end();
             it++) {
            extracted_promoters.push_back(promoters[it->second]);
            promoters.erase(it->second);
        }

        prom_pos.erase(first, prom_pos.end());

        // Find the last promoters in the interval
        auto end = prom_pos.lower_bound(pos_2);


        // Extract the promoters (remove them from the individual's list and put them in extracted_promoters)
        for (auto it = prom_pos.begin();
             it != end;
             it++) {
            extracted_promoters.push_back(promoters[it->second]);
            promoters.erase(it->second);
        }

        prom_pos.erase(prom_pos.begin(),end);

    } else {

        auto first = prom_pos.lower_bound(pos_1);

        if (first == prom_pos.end() or first->first >= pos_2) {
            return;
        }

        // Find the last promoters in the interval
        auto end = prom_pos.lower_bound(pos_2);


        // Extract the promoters (remove them from the individual's list and put them in extracted_promoters)
        for (auto it = first;
             it != end;
             it++) {
            extracted_promoters.push_back(promoters[it->second]);
            promoters.erase(it->second);
        }

        prom_pos.erase(first, end);
    }
}

void Organism::promoters_included_in(int32_t pos_1, int32_t pos_2, std::list<Promoter*>& promoters_list) {
    if (pos_1 < pos_2) {
        int32_t seg_length = pos_2 - pos_1;

        if (seg_length >= PROM_SIZE) {
            lst_promoters(BETWEEN, pos_1, pos_2 - PROM_SIZE + 1,
                          promoters_list);
        }
    }
    else {
        int32_t seg_length = dna_->length() + pos_2 - pos_1;
        if (seg_length >= PROM_SIZE) {
            bool is_near_end_of_genome = (pos_1 + PROM_SIZE > dna_->length());
            bool is_near_beginning_of_genome = (pos_2 - PROM_SIZE < 0);

            if (!is_near_end_of_genome && !is_near_beginning_of_genome) {
                lst_promoters(AFTER, pos_1, -1, promoters_list);
                lst_promoters(BEFORE, -1, pos_2 - PROM_SIZE + 1,
                              promoters_list);
            }
            else if (!is_near_end_of_genome) // => && is_near_beginning_of_genome
            {
                lst_promoters(BETWEEN, pos_1, pos_2 - PROM_SIZE + 1 +
                                                       dna_->length(),
                              promoters_list);
            }
            else if (!is_near_beginning_of_genome) // => && is_near_end_of_genome
            {
                lst_promoters(AFTER, pos_1, -1, promoters_list);
                lst_promoters(BEFORE, -1, pos_2 - PROM_SIZE + 1,
                              promoters_list);
            }
            else // is_near_end_of_genome && is_near_beginning_of_genome
            {
                lst_promoters(BETWEEN, pos_1, pos_2 - PROM_SIZE + 1 +
                                                       dna_->length(),
                              promoters_list);
            }
        }
    }
}


void Organism::lst_promoters(Position before_after_btw, // with regard to the strand's reading direction
                                         int32_t pos1,
                                         int32_t pos2,
                                         std::list<Promoter*>& promoters_list) {
    auto it_begin = prom_pos.begin();
    auto it_end = prom_pos.end();

    if (before_after_btw != BEFORE && pos1 != -1) {
            auto tmp_it = prom_pos.lower_bound(pos1);
            if (tmp_it == prom_pos.end())
                return;

            if (tmp_it!=prom_pos.end()) it_begin = tmp_it;

    }

    if (before_after_btw != AFTER && pos2 != -1) {
            auto tmp_it = prom_pos.lower_bound(pos2);

            if (tmp_it!=prom_pos.end()) it_end = tmp_it;
    }

    for (auto it = it_begin; it!=it_end; it++) {
        promoters_list.push_back(promoters[it->second]);
    }
}
