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




#ifndef AEVOL_JUMPING_MT_H_
#define AEVOL_JUMPING_MT_H_


// =================================================================
//                              Libraries
// =================================================================
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <zlib.h>

#include <vector>


// =================================================================
//                            Project Files
// =================================================================
#include "SFMT-src-1.4/SFMT.h"
#include "SFMT-src-1.4/jump/SFMT-jump.h"




// =================================================================
//                          Class declarations
// =================================================================




// MT_RAND_MAX = 2^32-1
#define MT_RAND_MAX         4294967295.0
#define MT_RAND_MAX_PLUS_1  4294967296.0

class JumpingMT
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    JumpingMT(const uint32_t& simple_seed);   // Initialize with a simple uint32_t
    JumpingMT(const JumpingMT &model);    // Create a copy of an existing generator
    JumpingMT(gzFile backup_file);           // Load from a gz backup file

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~JumpingMT();

    // =================================================================
    //                        Accessors: getters
    // =================================================================

    // =================================================================
    //                        Accessors: setters
    // =================================================================

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    inline double   random();         // Double in [0, 1[ (uniform distribution)
    inline int8_t   random(int8_t max);   // ~
    inline int16_t  random(int16_t max);  // ~
    inline int32_t  random(int32_t max);  // ~ > Integer in [0, max[ (uniform distribution)
    inline int64_t  random(int64_t max);  // ~
    int32_t         binomial_random(int32_t nb, double prob); // Binomial drawing of parameters (nb, prob)
    double          gaussian_random();                    // Double following a Standard Normal distribution
    int32_t          roulette_random(double* probs, int32_t nb_elts, bool verbose = false); // Roulette selection
    void            multinomial_drawing (int32_t* destination, double* source, int32_t nb_drawings, int32_t colors);
    // Multinomial drawing of parameters (nb, {source[0], source[1], ... source[colors-1]})

    void jump();

    void save(gzFile backup_file) const;

    // =================================================================
    //                           Public Attributes
    // =================================================================
    static int32_t nb_jumps;
    static double  jump_time;
/*



    std::vector<double> pickones;
    std::vector<double> pickones2;
    std::vector<double> pickones3;

    std::vector<double> cloned_probs;*/
  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    JumpingMT()
    {
      printf("%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    };
    /*JumpingMT(const JumpingMT &model)
    {
      printf("%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    };*/


    // =================================================================
    //                           Protected Methods
    // =================================================================
    static double gammln(double X);

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    sfmt_t* sfmt_;

};


// =====================================================================
//                           Getters' definitions
// =====================================================================

// =====================================================================
//                           Setters' definitions
// =====================================================================

// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================
/*!
  Draw a double precision real-number in [0, 1) with a uniform distribution
 */
inline double JumpingMT::random()
{
  return sfmt_genrand_real2(sfmt_);
}

/*!
  Draw an 8-bit integer in [0, max[ with a uniform distribution
 */
inline int8_t JumpingMT::random(int8_t max)
{
  return (int8_t) floor(((double)max) * sfmt_genrand_real2(sfmt_));
}

/*!
  Draw an 16-bit integer in [0, max[ with a uniform distribution
 */
inline int16_t JumpingMT::random(int16_t max)
{
  return (int16_t) floor(((double)max) * sfmt_genrand_real2(sfmt_));
}

/*!
  Draw an 32-bit integer in [0, max[ with a uniform distribution
 */
inline int32_t JumpingMT::random(int32_t max)
{
  return (int32_t) floor(((double)max) * sfmt_genrand_real2(sfmt_));
}

/*!
  Draw an 64-bit integer in [0, max[ with a uniform distribution
 */
inline int64_t JumpingMT::random(int64_t max)
{
  return (int64_t) floor(((double)max) * sfmt_genrand_real2(sfmt_));
}

#endif // AEVOL_JUMPING_MT_H_
