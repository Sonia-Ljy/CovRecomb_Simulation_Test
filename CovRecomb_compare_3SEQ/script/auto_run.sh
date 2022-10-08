#!/bin/bash
cd /home/soniali/Desktop/CovRecomb-Global-Version/Simulation_Test/Simulator_compare3SEQ_renew

echo "\n Start for simulation datasets generation and Comparison between the 3SEQ and CovRecomb method.\n"

# Ideal simulation, with less number of differential feature mutations; Figures S3A-B
for filepath in $*"/home/soniali/Desktop/03_CovRecomb/CovRecomb_Simulation_Test/CovRecomb_compare_3SEQ/"
do
	python Simulator_CovRecombTest.py -smp 6 -sr 1 -f filepath
	wait
	python seq3_run.py -smp 6 -sr 1 -f filepath
	wait
	python Compare_CoverageRate.py -smp 6 -sr 1 -f filepath
	wait
	echo "\n End of test. \n"
done

echo `date`

# Ideal simulation, with more number of differential feature mutations; Figures S3C-D
for filepath in $*"/home/soniali/Desktop/03_CovRecomb/CovRecomb_Simulation_Test/CovRecomb_compare_3SEQ/"
do
	python3 Simulator_CovRecombTest.py -smp 8 -sr 1 -f filepath
	wait
	python3 seq3_run.py -smp 8 -sr 1 -f filepath
	wait
	python3 Compare_CoverageRate.py -smp 8 -sr 1 -f filepath
	wait
	echo "\n End of test. \n"
done

echo `date`

# Analog simulation, with less number of differential feature mutations; Figures S3E-F
for filepath in $*"/home/soniali/Desktop/03_CovRecomb/CovRecomb_Simulation_Test/CovRecomb_compare_3SEQ/"
do
	python3 Simulator_CovRecombTest.py -smp 6 -sr 10 -f filepath
	wait
	python3 seq3_run.py -smp 6 -sr 10 -f filepath
	wait
	python3 Compare_CoverageRate.py -smp 6 -sr 10 -f filepath
	wait
	echo "\n End of test. \n"
done

echo `date`

# Analog simulation, with more number of differential feature mutations; Figures S3G-H
do
	python3 Simulator_CovRecombTest.py -smp 8 -sr 10 -f filepath
	wait
	python3 3seq_run.py -smp 8 -sr 10 -f filepath
	wait
	python3 Compare_CoverageRate.py -smp 8 -sr 10 -f filepath
	wait
	echo "\n End of test. \n"
done

echo `date`
