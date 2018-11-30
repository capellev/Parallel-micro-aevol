// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
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
//*****************************************************************************




// =================================================================
//                              Includes
// =================================================================
#include <inttypes.h>
#include <cstring>

#include <list>
#include <vector>
#include <memory>

#include <gtest/gtest.h>

#include "Organism.h"
#include "ExpManager.h"

using std::list;
using std::vector;




//############################################################################
//                                                                           #
//                         Class IndividualTest                              #
//                                                                           #
//############################################################################
class IndividualTest : public testing::Test {
 protected:
  virtual void SetUp(void);
  virtual void TearDown(void);

  std::shared_ptr<Organism> indiv1;
  std::shared_ptr<Organism> indiv2;
};

// ===========================================================================
//                                 Public Methods
// ===========================================================================
void IndividualTest::SetUp(void)
{
  auto w_max = 0.033333333;
  auto selection_pressure = 2000;
  ExpManager* exp_m = new ExpManager(1, 2, 1, 0.000001, 5000, w_max, selection_pressure, 1000);

  for (int i = 0; i <= 1; ++i) {
    // Build ad-hoc genomes
    // (and reverse to test the same things on the lagging strand.):
    //
    // indiv1: (AS + prom + AS + AG + AS + term + AS + prom + AS)
    // indiv2: (AS + AG + AS + term + AS + prom + AS)
    //
    // AS = Arbitrary Sequence
    // AG = Arbitrary Gene
    // Do not modify the sequences !

    // Define a few arbitrary sequences
    char as[5][10] = {
        "0011",
        "11101",
        "110011",
        "11000",
        "000101"
    };// len = 26

    // Define an arbitrary gene
    char gene[255];
    sprintf(gene, "0110110011000100110110010001");// len= 28 (gene = 4AA)

    // Define an arbitrary terminator
    char term[11 + 1] = "01000001101"; //len = 11

    // Define a few arbitrary promoters
    char prom[2][23] = {
        "0101010001110110010110", // dist from consensus: 2 => basal level: 0.6
        "0101011001110010010010"  // dist from consensus: 1 => basal level: 0.8
    }; // len = 22

    // Construct a gene with these arbitrary sequences
    auto genome_size = 1024;
    char* genome = new char[genome_size];
    sprintf(genome, "%s%s%s%s%s%s%s%s%s", as[0], prom[0], as[1], gene, as[2],
            term, as[3], prom[1], as[4]);
    // 0 AS0 4 PROM0 26 AS1 31 GENE 59 AS2 65 TERM 76 AS3 81 PROM1 103 AS4 109(==0)

    // Build indiv1
    int indiv_id = 0;

    exp_m->internal_organisms_[indiv_id] = std::make_shared<Organism>(exp_m, genome, indiv_id);
    indiv1 = exp_m->internal_organisms_[indiv_id];

    // Do transcription and translation
    exp_m->start_stop_RNA(indiv_id);
    exp_m->compute_RNA(indiv_id);

    exp_m->start_protein(indiv_id);
    exp_m->compute_protein(indiv_id);

    exp_m->translate_protein(indiv_id, w_max);

    exp_m->compute_phenotype(indiv_id);

    exp_m->compute_fitness(indiv_id, selection_pressure);



    // Build indiv2
    ++indiv_id;
    genome = new char[genome_size];
    sprintf(genome, "%s%s%s%s%s%s%s", as[0], gene, as[1], term, as[2], prom[1],
            as[3]);
    // 0 AS0 4 GENE 32 AS1 37 TERM 48 AS2 54 PROM1 76 AS3 82

    exp_m->internal_organisms_[indiv_id] = std::make_shared<Organism>(exp_m, genome, indiv_id);
    indiv2 = exp_m->internal_organisms_[indiv_id];

    // Do transcription and translation
    exp_m->start_stop_RNA(indiv_id);
    exp_m->compute_RNA(indiv_id);

    exp_m->start_protein(indiv_id);
    exp_m->compute_protein(indiv_id);

    exp_m->translate_protein(indiv_id, w_max);

    exp_m->compute_phenotype(indiv_id);

    exp_m->compute_fitness(indiv_id, selection_pressure);
  }


    // ***************************************************************************
    // The following commented code allows to print stuff about rnas and proteins

    // printf("%"PRId32" rnas and %"PRId32" prots\n", indiv4->rna_list()->size(), indiv4->protein_list()->size());
    // list_node<Rna*>* rna_node = indiv4->rna_list()->first();
    // while (rna_node != NULL)
    // {
    //   printf("%s rna at pos %"PRId32" (%f, %"PRId32")\n",
    //           rna_node->obj()->strand() == LEADING ? "LEADING":"LAGGING",
    //           rna_node->obj()->promoter_pos(),
    //           rna_node->obj()->basal_level(),
    //           rna_node->obj()->transcript_length());

    //   rna_node = rna_node->next();
    // }

    // list_node<Protein*>* protein_node = indiv4->protein_list()->first();
    // while (protein_node != NULL)
    // {
    //   printf("%s protein at pos %"PRId32" (length: %"PRId32", concentr: %f, nb_rnas: %"PRId32")\n",
    //           protein_node->obj()->strand() == LEADING ? "LEADING":"LAGGING",
    //           protein_node->obj()->shine_dal_pos(),
    //           protein_node->obj()->length(),
    //           protein_node->obj()->concentration(),
    //           protein_node->obj()->rna_list()->size());

    //   protein_node = protein_node->next();
    // }
}

void IndividualTest::TearDown(void)
{
}


TEST_F(IndividualTest, TestIndiv1)
{
  // Check that we have the right number of promoters, terminators etc
  // and at the right positions
  // "right" means those values we have computed by hand

  // Check genome size
  EXPECT_EQ(109, indiv1->length());

  // Check RNA list
  const auto& rna_list = indiv1->rnas;
  EXPECT_EQ(2, rna_list.size());
  auto rna = rna_list.front();
  EXPECT_EQ(4, rna->begin);
  EXPECT_FLOAT_EQ(0.6, rna->e);
  EXPECT_EQ(50, rna->length);

  rna = rna_list.back();
  EXPECT_EQ(81, rna->begin);
  EXPECT_FLOAT_EQ(0.8, rna->e);
  EXPECT_EQ(82, rna->length);



  // Check protein list
  auto prot_list = indiv1->proteins;
  EXPECT_EQ(2, prot_list.size());
  auto prot = prot_list.front();
  EXPECT_EQ(31, prot->protein_start);
  EXPECT_EQ(4, prot->protein_length);
  EXPECT_FLOAT_EQ(1.4, prot->e);
  prot = prot_list.back();
  EXPECT_EQ(31, prot->protein_start);
  EXPECT_EQ(4, prot->protein_length);
  EXPECT_FLOAT_EQ(0.8, prot->e);
  EXPECT_FALSE(prot->is_init_);
}

TEST_F(IndividualTest, TestIndiv2)
{
  // Check genome size
  EXPECT_EQ(81, indiv2->length());

  // Check RNA list
  auto rna_list = indiv2->rnas;
  EXPECT_EQ(1, rna_list.size());
  auto rna = rna_list.front();
  EXPECT_EQ(54, rna->begin);
  EXPECT_FLOAT_EQ(0.8, rna->e);
  EXPECT_EQ(42, rna->length);

  // Check protein list
  const auto& prot_list = indiv2->proteins;
  EXPECT_EQ(1, prot_list.size());
  Protein* prot = prot_list.front();
  EXPECT_EQ(4, prot->protein_start);
  EXPECT_EQ(4, prot->protein_length);
  EXPECT_FLOAT_EQ(0.8, prot->e);
}

// ===========================================================================
//                                Protected Methods
// ===========================================================================

// ===========================================================================
//                              Non inline accessors
// ===========================================================================
