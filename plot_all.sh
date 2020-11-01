#!/bin/bash
export DYLD_LIBRARY_PATH=/usr/local/boost_1_67_0/stage/lib
#cd /Users/Kianusch/Documents/Studium/Semester/WiSe1819/Bachelor-Arbeit/Heavy-Ion-Collision-Analysis

#Fig. 1
python3 plot_hist.py

#Fig. 2
python3 plot_W_comparison.py

#Fig. 3
python3 plot_one_point_centrality.py

#Fig. 4
python3 plot_one_point_comparison_centrality.py

#python3 plot_two_point_comparison_OnePoint.py

#Fig. 5
python3 plot_two_point_comparison_paper_0-1_horizontal.py

#Fig. 6
python3 plot_two_point_comparison_paper_20-21_horizontal.py

#Fig. 7
python3 plot_two_point_IPSM.py

#Fig. 8
python3 plot_two_point_comparison_centrality.py

#python3 plot_two_point_comparison_universality.py





