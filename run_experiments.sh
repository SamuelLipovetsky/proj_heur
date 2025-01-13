
echo 'file;alfa_param;cutoff_time;greedy_mege;nuqcl;random_seed'> results.txt
for file in ./data/*
do
 
    for alfa_param in 0.91 0.95 0.98
    do
        
        for cutoff_time in 3
        do
            random_seed=$RANDOM 
            echo -n "$file;$alfa_param;$cutoff_time;"
            ./proj_heuristic/quasi_clique_finder $file $cutoff_time $alfa_param $random_seed
            echo -n ";"
            ./base_heuristic/nuqclq $file $cutoff_time $alfa_param $random_seed
            echo -n ";$random_seed"
            echo ";"

        done
    done 
done >> results.txt