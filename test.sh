make -C build
rm -f max_clique && ln -s build/max_clique max_clique
time ./max_clique < anna.col
time ./max_clique < brock200_2.clq
time ./max_clique < brock200_4.clq
time ./max_clique < p_hat300-2.clq
