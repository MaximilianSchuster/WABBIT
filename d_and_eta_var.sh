# change into d_pip_directory, execute wabbit
cd d_pip_0_2
mpirun -n 4 ./wabbit plenum.ini --memory=4gb
zip d_pip_0_2_080219.zip *.h5
echo "1 out of 5 done"
cd ..
cd d_pip_0_1
mpirun -n 4 ./wabbit plenum.ini --memory=4gb
zip d_pip_0_1_080219.zip *.h5
echo "2 out of 5 done"

cd ..
cd d_pip_0_05
mpirun -n 4 ./wabbit plenum.ini --memory=4gb
zip d_pip_0_05_080219.zip *.h5
echo "3 out of 5 done"

cd ..
cd eta_e6
mpirun -n 4 ./wabbit plenum.ini --memory=4gb
zip eta_e6_080219.zip *.h5
echo "4 out of 5 done"

cd ..
cd eta_e5
mpirun -n 4 ./wabbit plenum.ini --memory=4gb
zip eta_e5_080219.zip *.h5
echo "5 out of 5 done"

