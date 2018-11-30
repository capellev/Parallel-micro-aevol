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



#ifndef RAEVOL_CUDA_MUTATIONEVENT_H
#define RAEVOL_CUDA_MUTATIONEVENT_H

#include <cstdint>
#include <vector>

enum MutationEventType {
    DO_SWITCH           = 0,
  NONE
};

/**
 * Mutation event class
 */
class Dna;

class MutationEvent {

 public:
    MutationEvent() = default;
    ~MutationEvent();

    void switch_pos(int32_t pos);

    int32_t type() { return type_; };

    int32_t pos_1() { return pos_1_; }
    int32_t pos_2() { return pos_2_; }
    int32_t pos_3() { return pos_3_; }

    int32_t number() { return number_; }

    int32_t invert() { return invert_; }


    Dna* seq() { return seq_; }


 private:
    int32_t type_;

    int32_t pos_1_,pos_2_,pos_3_;

    int32_t number_; // insertion or deletion

    bool invert_;


    Dna* seq_;

};



#endif //RAEVOL_CUDA_MUTATIONEVENT_H
