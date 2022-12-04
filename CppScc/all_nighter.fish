function all_nighter --argument times --description="A script to take build all implementations and run them for different thread counts, expected to take many hours"
    cd ~/Desktop/pdsProject/CppScc

    mkdir -p results/serial results/OpenMP results/OpenCilk results/pthread

    # Make measurements for the serial version, only one measurement needed
    echo "Serial"
    cd serial
    make clean
    make
    ./colorSCC ../../matrices $times 0 1 | tee ../results/serial/serial_cold.txt
    cd ..

    # Make measurements for the OpenMP version
    echo "OpenMP"
    cd OpenMP
    make clean
    make
    for i in (seq 2 2 12)
        export OMP_NUM_THREADS=$i
        echo "Running with $i threads"
        ./colorSCC ../../matrices $times 0 1 | tee ../results/OpenMP/threads_$i.txt
    end
    cd ..

    # run the serial again to observe thermal throttling effects
    echo "Serial"
    cd serial
    make clean
    make
    ./colorSCC ../../matrices $times 0 1 | tee ../results/serial/serial_hot.txt
    cd ~/Desktop/pdsProject/CppScc
 
    # Make measurements for the OpenCilk version
    echo "OpenCilk"
    cd OpenCilk
    make clean
    make
    for i in (seq 2 2 12)
        export CILK_NWORKERS=$i
        echo "Running with $i threads"
        ./colorSCC ../../matrices $times 0 1 | tee ../results/OpenCilk/threads_$i.txt
    end
    cd ..

    # Make measurements for the pthread version
    echo "Pthread"
    cd pthread
    make clean
    make
    for i in (seq 2 2 12)
        echo "Running with $i threads"
        ./colorSCC ../../matrices $times 0 1 $i | tee ../results/pthread/threads_$i.txt
    end
    cd ..



    #IF done before i check the results, keep it up for more threads

    # Make measurements for the OpenMP version
    echo "OpenMP"
    cd OpenMP
    make clean
    make
    for i in (seq 12 2 24)
        export OMP_NUM_THREADS=$i
        echo "Running with $i threads"
        ./colorSCC ../../matrices $times 0 1 | tee ../results/OpenMP/threads_$i.txt
    end
    cd ..

    # Make measurements for the OpenCilk version
    echo "OpenCilk"
    cd OpenCilk
    make clean
    make
    for i in (seq 12 2 24)
        export CILK_NWORKERS=$i
        echo "Running with $i threads"
        ./colorSCC ../../matrices $times 0 1 | tee ../results/OpenCilk/threads_$i.txt
    end
    cd ..

    # Make measurements for the pthread version
    echo "Pthread"
    cd pthread
    make clean
    make
    for i in (seq 12 2 24)
        echo "Running with $i threads"
        ./colorSCC ../../matrices $times 0 1 $i | tee ../results/pthread/threads_$i.txt
    end
    cd ..
    

   #ring the alarm and wake me up please
end