# 全局规划器的一些说明

## 算法概述

* 总体算法过程可描述如下

输入：代价地图costmap，起点start，终点goal

输出：path，数据结构为PoseStamped array

1. Dji expander计算potential

1. gradient path_maker从potential中得到路径

## 1. Dji Expander

输入：代价地图，起点，终点

输出：potential（至今还不太清楚potential是什么意思，可能是Dji算法中的优先级吧）

### 1.1 wavefront algorithm

* wavefront 算法本来是人工势场法中的一种，现在作为
