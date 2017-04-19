/*********************************************************************
 *
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2008, 2013, Willow Garage, Inc.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Willow Garage, Inc. nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * Author: Eitan Marder-Eppstein
 *         David V. Lu!!
 *********************************************************************/
#ifndef _DIJKSTRA_H
#define _DIJKSTRA_H

#define PRIORITYBUFSIZE 10000
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include <global_planner/planner_core.h>
#include <global_planner/expander.h>

// inserting onto the priority blocks
#define push_cur(n)  { if (n>=0 && n<ns_ && !pending_[n] && getCost(costs, n)<lethal_cost_ && currentEnd_<PRIORITYBUFSIZE){ currentBuffer_[currentEnd_++]=n; pending_[n]=true; }}
#define push_next(n) { if (n>=0 && n<ns_ && !pending_[n] && getCost(costs, n)<lethal_cost_ &&    nextEnd_<PRIORITYBUFSIZE){    nextBuffer_[   nextEnd_++]=n; pending_[n]=true; }}
#define push_over(n) { if (n>=0 && n<ns_ && !pending_[n] && getCost(costs, n)<lethal_cost_ &&    overEnd_<PRIORITYBUFSIZE){    overBuffer_[   overEnd_++]=n; pending_[n]=true; }}

namespace global_planner {
class DijkstraExpansion : public Expander {
    public:
        DijkstraExpansion(PotentialCalculator* p_calc, int nx, int ny);
        ~DijkstraExpansion();
        // 计算potentials，根据起点和终点，以及代价值，结合加入的neutral_cost等参数，进行wavefront处理，然后使用dji方法计算路径，返回potential指针
        bool calculatePotentials(unsigned char* costs, double start_x, double start_y, double end_x, double end_y, int cycles,
                                float* potential);

        /**
         * @brief  Sets or resets the size of the map
         * @param nx The x size of the map
         * @param ny The y size of the map
         */
        void setSize(int nx, int ny); /**< sets or resets the size of the map */

        void setNeutralCost(unsigned char neutral_cost) {
            neutral_cost_ = neutral_cost;
            priorityIncrement_ = 2 * neutral_cost_;
        }

        void setPreciseStart(bool precise){ precise_ = precise; }
    private:

        /**
         * @brief  Updates the cell at index n
         * @param costs The costmap
         * @param potential The potential array in which we are calculating
         * @param n The index to update
         */
        // 更新栅格？为什么要更新？输入为代价，输出为potential
        void updateCell(unsigned char* costs, float* potential, int n); /** updates the cell at index n */

        // 获取代价值，最高为lethal_cost，即只对区间内部的代价值进行处理c = c * factor_ + neutral_cost_;
        // 处理完毕后还进行修整，如果超出值域则
        float getCost(unsigned char* costs, int n) {
            float c = costs[n];
            // 如果
            if (c < lethal_cost_ - 1 || (unknown_ && c==255)) {
                c = c * factor_ + neutral_cost_;
                if (c >= lethal_cost_)
                    c = lethal_cost_ - 1;
                return c;
            }
            return lethal_cost_;
        }

        /** block priority buffers */
        int *buffer1_, *buffer2_, *buffer3_; /**< storage buffers for priority blocks */
        int *currentBuffer_, *nextBuffer_, *overBuffer_; /**< priority buffer block ptrs */
        int currentEnd_, nextEnd_, overEnd_; /**< end points of arrays */
        bool *pending_; /**< pending_ cells during propagation */
        bool precise_;

        /** block priority thresholds */
        float threshold_; /**< current threshold */
        float priorityIncrement_; /**< priority threshold increment */

};
} //end namespace global_planner
#endif
