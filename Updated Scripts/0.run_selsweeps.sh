echo "run slim simulations for SNM and for selective sweeps: Complete selective sweep, incomplete selective sweep, soft selective sweep:"
echo "conditions of the SLiM simulations are included in the script sweeps_ind.slim: mu=2e-7, r=1e-7, L=5e5, 2Ne=1000"
rm *.ms
rm *.ms.log
NITER=100
for i in $(seq 1 ${NITER})
do
    echo iteration $i/$NITER
    slim -t -m -d "freq_sel_init=0.001" -d "freq_sel_end=0.999" -d "s_beneficial=0.00" -d "ind_sample_size=50" -d "file_output1='./slim.snm.output_all.ms'" ./sweeps_ind.slim 1>>  ./slim.snm.output_all.ms.log
    slim -t -m -d "freq_sel_init=0.001" -d "freq_sel_end=1.000" -d "s_beneficial=0.05" -d "ind_sample_size=50" -d "file_output1='./slim.selsweep.output_all.ms'" ./sweeps_ind.slim 1>>  ./slim.selsweep.output_all.ms.log
    slim -t -m -d "freq_sel_init=0.001" -d "freq_sel_end=0.750" -d "s_beneficial=0.05" -d "ind_sample_size=50" -d "file_output1='./slim.incomplete_selsweep.output_all.ms'" ./sweeps_ind.slim 1>> ./slim.incomplete_selsweep.output_all.ms.log
    slim -t -m -d "freq_sel_init=0.500" -d "freq_sel_end=1.000" -d "s_beneficial=0.05" -d "ind_sample_size=50" -d "file_output1='./slim.standing_selsweep.output_all.ms'" ./sweeps_ind.slim 1>>  ./slim.standing_selsweep.output_all.ms.log
done
#rm tmp_slim_*.txt

// Keywords: Selevted sweep at 5e5 position, Complete and incomplete sweep

//    slim -t -m -d "freq_sel_init=0.5" -d "freq_sel_end=0.6" -d "s_beneficial=0.0001" -d "ind_sample_size=20" -d "file_output1='./slim.snm.output_file.ms'" ./sweeps_ind.slim

//   slim -t -m -d "freq_sel_init=0.05" -d "freq_sel_end=1.0" -d "s_beneficial=0.05" -d "ind_sample_size=20" -d "file_output1='./slim.selsweep.output_file.ms'" ./sweeps_ind.slim

//    slim -t -m -d "freq_sel_init=0.05" -d "freq_sel_end=0.5" -d "s_beneficial=0.05" -d "ind_sample_size=20" -d "file_output1='./slim.incomplete_selsweep.output_file.ms'" ./sweeps_ind.slim

//    slim -t -m -d "freq_sel_init=0.5" -d "freq_sel_end=1.0" -d "s_beneficial=0.05" -d "ind_sample_size=20" -d "file_output1='./slim.standing_selsweep.output_file.ms'" ./sweeps_ind.slim

initialize() {
	initializeMutationRate(2e-7);
	initializeMutationType("m1", 0.5, "f", 0.0);
	initializeGenomicElementType("g1", m1, 1.0);
	initializeGenomicElement(g1, 0, 499999);
	initializeRecombinationRate(1e-7);
	
	defineConstant("Ne",500);
	
	if (exists("slimgui")) {
		defineConstant("freq_sel_init",0.001);
		defineConstant("freq_sel_end",0.999);
		defineConstant("s_beneficial",0.05);
		defineConstant("ind_sample_size",50);
		defineConstant("file_output1","~/slim.test_selsweep.output_file.ms");
	}
}
//do initialization of the population during 10000 generations (20*Ne) to have a population in equlibrium mutation-drift
1 early() { 
	// save this run's identifier, used to save and restore
	defineConstant("simID", getSeed());
	sim.addSubpop("p1", Ne); 
}
//introduce selection at a given posotopn using some conditions (frequency)
10000 late() {
	if(s_beneficial != 0) {	
		// save the state of the simulation
		sim.outputFull("./tmp_slim_" + simID + ".txt");	
		//target = sample(p1.genomes, 1);
		//target.addNewDrawnMutation(m2, 500000);
		muts = sim.mutations;
		//mut = muts[(sim.mutationFrequencies(p1, muts) == freq_sel_init)];
		mut = muts[(sim.mutationFrequencies(p1, muts) >= freq_sel_init) & (sim.mutationFrequencies(p1, muts) <= freq_sel_init + (freq_sel_end - freq_sel_init)/50)];
		
		if (size(mut))
		{
			mut = sample(mut, 1); //collect a single mutation
			mut.setSelectionCoeff(s_beneficial); //add beneficial slection
			pos_mut = mut.id; //look at the variant in the total variants
			print("position:" + mut.position);
		}
		else
		{
			cat("No contender of sufficient frequency found.\n");
		}
	}
	else {
			cat("ACHIEVED\n");
			allIndividuals = sim.subpopulations.individuals;
			//collect genomes from the SAME INDIVIDUAL
			sampledIndividuals = sample(allIndividuals, ind_sample_size,replace=F);
			g_1 = sampledIndividuals.genomes;
			m = sortBy(unique(g_1.mutations),"position");
			// print the number of segregating sites
			writeFile(filePath=file_output1,contents=("//"),append=T);
			writeFile(filePath=file_output1,contents=("segsites: "+ size(m)),append=T);
			//print the positions
			positions = format("%.6f", m.position / sim.chromosome.lastPosition);
			writeFile(filePath=file_output1,contents=("positions: "+ paste(positions, sep=" ")),append=T);
			//print the sampled genomes
			for (genome in g_1){
				hasMuts = (match(m,genome.mutations) >= 0);
				writeFile(filePath=file_output1,contents=(paste(asInteger(hasMuts),sep="")),append=T);
			}
			//sim.outputFull();
			print("Simulation finished");
			//finish simulation
			system("rm ./tmp_slim_" + simID + ".txt");
			sim.simulationFinished();
	}
}
//finish simulation if selective position achieve final frequency, otherwise repeat process from equilibrium saved population
10001:15000 late() {	
	if (sum(sim.mutations.selectionCoeff) == 0)
	{
		if (sum(sim.substitutions.selectionCoeff) == 0.0) {
			cat("Sweep mutation lost in gen. " + sim.cycle + "\n");
			// go back to initial saved generation 10000
			sim.readFromPopulationFile("./tmp_slim_" + simID + ".txt");			
			// start a newly seeded run
			setSeed(rdunif(1, 0, asInteger(2^62) - 1));
			muts = sim.mutations;
			//mut = muts[(sim.mutationFrequencies(p1, muts) == freq_sel_init)];
			mut = muts[sim.mutationFrequencies(p1, muts) > freq_sel_init & (sim.mutationFrequencies(p1, muts) <= freq_sel_init + (freq_sel_end - freq_sel_init)/50)];
			
			if (size(mut))
			{
				 mut = sample(mut, 1); //collect a single mutation
				 mut.setSelectionCoeff(s_beneficial); //add beneficial slection
				 pos_mut = mut.id; //look at the variant in the total variants
				 print("position:" + mut.position);
			}
			else
			{
				cat("No contender of sufficient frequency found.\n");
				//sim.outputFull();
				print("Simulation finished");
				//finish simulation
				system("rm ./tmp_slim_" + simID + ".txt");
				sim.simulationFinished();
			}
		}
		else {
			cat("Sweep mutation reached fixation.\n");
			allIndividuals = sim.subpopulations.individuals;
			//collect genomes from the SAME INDIVIDUAL
			sampledIndividuals = sample(allIndividuals, ind_sample_size,replace=F);
			g_1 = sampledIndividuals.genomes;
			m = sortBy(unique(g_1.mutations),"position");
			// print the number of segregating sites
			writeFile(filePath=file_output1,contents=("//"),append=T);
			writeFile(filePath=file_output1,contents=("segsites: "+ size(m)),append=T);
			//print the positions
			positions = format("%.6f", m.position / sim.chromosome.lastPosition);
			writeFile(filePath=file_output1,contents=("positions: "+ paste(positions, sep=" ")),append=T);
			//print the sampled genomes
			for (genome in g_1){
				hasMuts = (match(m,genome.mutations) >= 0);
				writeFile(filePath=file_output1,contents=(paste(asInteger(hasMuts),sep="")),append=T);
			}
			//sim.outputFull();
			print("Simulation finished");
			//finish simulation
			system("rm ./tmp_slim_" + simID + ".txt");
			sim.simulationFinished();
		}
	}
	//monitorize the beneficial mutation
	muts = sim.mutations;
	pos_mut = which(sim.mutations.selectionCoeff > 0);
	if(size(pos_mut)) {
		print(sim.mutationFrequencies(p1,muts[pos_mut])); 
		if (sim.mutationFrequencies(p1,muts[pos_mut]) >= freq_sel_end) {	
			print("freq_end = " + sim.mutationFrequencies(p1,muts[pos_mut]));
			cat("ACHIEVED\n");
			allIndividuals = sim.subpopulations.individuals;
			//collect genomes from the SAME INDIVIDUAL
			sampledIndividuals = sample(allIndividuals, ind_sample_size,replace=F);
			g_1 = sampledIndividuals.genomes;
			m = sortBy(unique(g_1.mutations),"position");
			// print the number of segregating sites
			writeFile(filePath=file_output1,contents=("//"),append=T);
			writeFile(filePath=file_output1,contents=("segsites: "+ size(m)),append=T);
			//print the positions
			positions = format("%.6f", m.position / sim.chromosome.lastPosition);
			writeFile(filePath=file_output1,contents=("positions: "+ paste(positions, sep=" ")),append=T);
			//print the sampled genomes
			for (genome in g_1){
				hasMuts = (match(m,genome.mutations) >= 0);
				writeFile(filePath=file_output1,contents=(paste(asInteger(hasMuts),sep="")),append=T);
			}
			//sim.outputFull();
			print("Simulation finished");
			//finish simulation
			system("rm ./tmp_slim_" + simID + ".txt");
			sim.simulationFinished();
		}	
	}
}

----------------------------------------------------------------------------------------------------------------------------------
echo "run slim simulations for SNM and for selective sweeps: Complete selective sweep, incomplete selective sweep, soft selective sweep:"
echo "conditions of the SLiM simulations are included in the script sweeps_ind.slim: mu=2e-6, r=1e-8, L=5e5, 2Ne=1000"
rm *.ms
rm *.ms.log
NITER=10
for i in $(seq 1 ${NITER})
do
    echo iteration $i/$NITER
    slim -t -m -d "freq_sel_init=0.001" -d "freq_sel_end=0.999" -d "s_beneficial=0.00" -d "ind_sample_size=50" -d "file_output1='./slim.snm.output_all.ms'" ./sweeps_ind.slim 1>>  ./slim.snm.output_all.ms.log
    slim -t -m -d "freq_sel_init=0.001" -d "freq_sel_end=1.000" -d "s_beneficial=0.5" -d "ind_sample_size=50" -d "file_output1='./slim.selsweep.output_all.ms'" ./sweeps_ind.slim 1>>  ./slim.selsweep.output_all.ms.log
    slim -t -m -d "freq_sel_init=0.001" -d "freq_sel_end=0.750" -d "s_beneficial=0.5" -d "ind_sample_size=50" -d "file_output1='./slim.incomplete_selsweep.output_all.ms'" ./sweeps_ind.slim 1>> ./slim.incomplete_selsweep.output_all.ms.log
    slim -t -m -d "freq_sel_init=0.500" -d "freq_sel_end=1.000" -d "s_beneficial=0.5" -d "ind_sample_size=50" -d "file_output1='./slim.standing_selsweep.output_all.ms'" ./sweeps_ind.slim 1>>  ./slim.standing_selsweep.output_all.ms.log
done
#rm tmp_slim_*.txt

-------------------------------------------------------------------------------------------------------------------------------------
echo "run slim simulations for SNM and for selective sweeps: Complete selective sweep, incomplete selective sweep, soft selective sweep:"
echo "conditions of the SLiM simulations are included in the script sweeps_ind.slim: mu=2e-7, r=1e-7, L=5e5, 2Ne=1000"
rm *.ms
rm *.ms.log
NITER=100
for i in $(seq 1 ${NITER})
do
    echo iteration $i/$NITER
    slim -t -m -d "freq_sel_init=0.001" -d "freq_sel_end=0.999" -d "s_beneficial=0.00" -d "ind_sample_size=50" -d "file_output1='./slim.snm.output_all.ms'" ./sweeps_ind.slim 1>>  ./slim.snm.output_all.ms.log
    slim -t -m -d "freq_sel_init=0.001" -d "freq_sel_end=1.000" -d "s_beneficial=0.50" -d "ind_sample_size=50" -d "file_output1='./slim.selsweep.output_all.ms'" ./sweeps_ind.slim 1>>  ./slim.selsweep.output_all.ms.log
    slim -t -m -d "freq_sel_init=0.001" -d "freq_sel_end=0.750" -d "s_beneficial=0.50" -d "ind_sample_size=50" -d "file_output1='./slim.incomplete_selsweep.output_all.ms'" ./sweeps_ind.slim 1>> ./slim.incomplete_selsweep.output_all.ms.log
    slim -t -m -d "freq_sel_init=0.500" -d "freq_sel_end=1.000" -d "s_beneficial=0.50" -d "ind_sample_size=50" -d "file_output1='./slim.standing_selsweep.output_all.ms'" ./sweeps_ind.slim 1>>  ./slim.standing_selsweep.output_all.ms.log
done
#rm tmp_slim_*.txt