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
#include<global_planner/dijkstra.h>
#include <algorithm>
namespace global_planner {

DijkstraExpansion::DijkstraExpansion(PotentialCalculator* p_calc, int nx, int ny) :
        Expander(p_calc, nx, ny), pending_(NULL), precise_(false) {
    // priority buffers
    // 哪里是优先级队列？就是普通数组呗
    buffer1_ = new int[PRIORITYBUFSIZE];
    buffer2_ = new int[PRIORITYBUFSIZE];
    buffer3_ = new int[PRIORITYBUFSIZE];

    priorityIncrement_ = 2 * neutral_cost_;
}

DijkstraExpansion::~DijkstraExpansion() {
  delete[] buffer1_;
  delete[] buffer2_;
  delete[] buffer3_;
  if (pending_)
      delete[] pending_;
}

//
// Set/Reset map size
//
void DijkstraExpansion::setSize(int xs, int ys) {
    Expander::setSize(xs, ys);
    if (pending_)
        delete[] pending_;

    pending_ = new bool[ns_];
    memset(pending_, 0, ns_ * sizeof(bool));
}

//
// main propagation function
// Dijkstra method, breadth-first
// runs for a specified number of cycles,
//   or until it runs out of cells to update,
//   or until the Start cell is found (atStart = true)
//
// dij方法，广度优先搜索，终止条件有三个：搜索超过指定数目循环次数，或者更新队列空，或者起始栅格(起点)被找到
// wavefront algorithm，给定代价值，起点终点，输出potential，这个potentian给path_maker沿着栅格或者梯度下降方法找到一个路径

bool DijkstraExpansion::calculatePotentials(unsigned char* costs, double start_x, double start_y, double end_x, double end_y,
                                           int cycles, float* potential) {
    cells_visited_ = 0;
    // priority buffers
    threshold_ = lethal_cost_;
    currentBuffer_ = buffer1_;
    currentEnd_ = 0;
    nextBuffer_ = buffer2_;
    nextEnd_ = 0;
    overBuffer_ = buffer3_;
    overEnd_ = 0;
    memset(pending_, 0, ns_ * sizeof(bool));
    std::fill(potential, potential + ns_, POT_HIGH);

    // set goal
    // 从栅格的x,y坐标映射到数组的index序号
    // 起点栅格是wavefront的终点
    int k = toIndex(start_x, start_y);

    // 是否为old_navfn_behavior,否则是precise=trut，起点为精确值
    // 老navfn行为是不精确的，新的方法考虑了精确的方面
    if(precise_)
    {
        // 八邻接方式
        // 考虑起点并不在栅格上，而是在栅格的边界上
        double dx = start_x - (int)start_x, dy = start_y - (int)start_y;
        dx = floorf(dx * 100 + 0.5) / 100;
        dy = floorf(dy * 100 + 0.5) / 100;
        potential[k] = neutral_cost_ * 2 * dx * dy;
        potential[k+1] = neutral_cost_ * 2 * (1-dx)*dy;
        potential[k+nx_] = neutral_cost_*2*dx*(1-dy);
        potential[k+nx_+1] = neutral_cost_*2*(1-dx)*(1-dy);//*/

        push_cur(k+2);
        push_cur(k-1);
        push_cur(k+nx_-1);
        push_cur(k+nx_+2);

        push_cur(k-nx_);
        push_cur(k-nx_+1);
        push_cur(k+nx_*2);
        push_cur(k+nx_*2+1);
    }else{
        // 将栅格的带价值压入优先级队列
        // 四邻接方式
        // 目标点处势场设置为0
        // #define push_cur(n)  { if (n>=0 && n<ns_ && !pending_[n] && getCost(costs, n)<lethal_cost_ && currentEnd_<PRIORITYBUFSIZE){ currentBuffer_[currentEnd_++]=n; pending_[n]=true; }}
        // #define push_next(n) { if (n>=0 && n<ns_ && !pending_[n] && getCost(costs, n)<lethal_cost_ &&    nextEnd_<PRIORITYBUFSIZE){    nextBuffer_[   nextEnd_++]=n; pending_[n]=true; }}
        // #define push_over(n) { if (n>=0 && n<ns_ && !pending_[n] && getCost(costs, n)<lethal_cost_ &&    overEnd_<PRIORITYBUFSIZE){    overBuffer_[   overEnd_++]=n; pending_[n]=true; }}
        // 三个优先级队列，每个队列有一个end，优先级为getcost的返回值，即经过加权处理的代价值
        // 目标点处势场设置为0        
        potential[k] = 0;
        push_cur(k+1); //右边邻接栅格
        push_cur(k-1); //左边邻接栅格
        push_cur(k-nx_); //上面邻接栅格
        push_cur(k+nx_); //下面邻接栅格
    }

    int nwv = 0;            // max priority block size
    int nc = 0;            // number of cells put into priority blocks
    int cycle = 0;        // which cycle we're on

    // set up start cell
    // 从这个栅格，即目标点栅格开始wavefront
    int startCell = toIndex(end_x, end_y);

    // for语句仅作为终止条件
    for (; cycle < cycles; cycle++) // go for this many cycles, unless interrupted
            {
        // 
        if (currentEnd_ == 0 && nextEnd_ == 0) // priority blocks empty
            return false;

        // stats
        nc += currentEnd_;
        if (currentEnd_ > nwv)
            nwv = currentEnd_;

        // reset pending_ flags on current priority buffer
        int *pb = currentBuffer_;
        int i = currentEnd_;
        while (i-- > 0)
            pending_[*(pb++)] = false;

        // process current priority buffer
        pb = currentBuffer_;
        i = currentEnd_;
        // 重要：对优先队列中的栅格进行更新
        while (i-- > 0)
            updateCell(costs, potential, *pb++);

        // swap priority blocks currentBuffer_ <=> nextBuffer_
        currentEnd_ = nextEnd_;
        nextEnd_ = 0;
        pb = currentBuffer_;        // swap buffers
        currentBuffer_ = nextBuffer_;
        nextBuffer_ = pb;

        // see if we're done with this priority level
        if (currentEnd_ == 0) {
            threshold_ += priorityIncrement_;    // increment priority threshold
            currentEnd_ = overEnd_;    // set current to overflow block
            overEnd_ = 0;
            pb = currentBuffer_;        // swap buffers
            currentBuffer_ = overBuffer_;
            overBuffer_ = pb;
        }

        // check if we've hit the Start cell
        if (potential[startCell] < POT_HIGH)
            break;
    }
    //ROS_INFO("CYCLES %d/%d ", cycle, cycles);
    if (cycle < cycles)
        return true; // finished up here
    else
        return false;
}

//
// Critical function: calculate updated potential value of a cell,
//   given its neighbors' values
// Planar-wave update calculation from two lowest neighbors in a 4-grid
// Quadratic approximation to the interpolated value
// No checking of bounds here, this function should be fast
//
// wavefront algorithms
// 输入为某一个栅格的经过加权处理的代价值，输出为一个经过更新的potential, potential更新是根据栅格临接栅格值来的。

// sqrt(2)/2
#define INVSQRT2 0.707106781
// sqrt(2)/2

inline void DijkstraExpansion::updateCell(unsigned char* costs, float* potential, int n) {
    cells_visited_++;

    // do planar wave update
    float c = getCost(costs, n);
    if (c >= lethal_cost_)    // don't propagate into obstacles
        return;
    
    // 这个函数，计算势场值。
    float pot = p_calc_->calculatePotential(potential, c, n);

    // now add affected neighbors to priority blocks
    if (pot < potential[n]) {
        float le = INVSQRT2 * (float)getCost(costs, n - 1);
        float re = INVSQRT2 * (float)getCost(costs, n + 1);
        float ue = INVSQRT2 * (float)getCost(costs, n - nx_);
        float de = INVSQRT2 * (float)getCost(costs, n + nx_);
        potential[n] = pot;
        //ROS_INFO("UPDATE %d %d %d %f", n, n%nx, n/nx, potential[n]);
        if (pot < threshold_)    // low-cost buffer block
                {
            if (potential[n - 1] > pot + le)
                push_next(n-1);
            if (potential[n + 1] > pot + re)
                push_next(n+1);
            if (potential[n - nx_] > pot + ue)
                push_next(n-nx_);
            if (potential[n + nx_] > pot + de)
                push_next(n+nx_);
        } else            // overflow block
        {
            if (potential[n - 1] > pot + le)
                push_over(n-1);
            if (potential[n + 1] > pot + re)
                push_over(n+1);
            if (potential[n - nx_] > pot + ue)
                push_over(n-nx_);
            if (potential[n + nx_] > pot + de)
                push_over(n+nx_);
        }
    }
}

} //end namespace global_planner
