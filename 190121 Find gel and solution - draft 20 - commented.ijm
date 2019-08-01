nBins = 200; //Number of bins in final ratio histogram
histMin = 0; //Minimum bin value for all histograms
histMax= 2; //Maximum bin value for all histograms
median = 1; //Size of 3D median filter to denoise raw stacks

// This version of the script uses the second derivative for all transition locating.
// It also reports the ratio of gel:soln above coverslip as well as ratio of gel:soln above gel


// Make sure results and log windows start fresh
if(isOpen("Results")){
	selectWindow("Results");
	run("Close");
}
if(isOpen("Log")){
	selectWindow("Log");
	run("Close");
}

// Batchmode = TRUE does not show images as analysis is proceeding - much faster
setBatchMode(true);

// Open directory and make a new output directory to store analysis
run("Bio-Formats Macro Extensions");
run("Close All");
dir = getDirectory("Choose input directory");
list = getFileList(dir);
outDir = dir + "Output_draft20/";
File.makeDirectory(outDir);
counter = 0;
newImage("Normalized Ratio Histogram", "32-bit black", nBins, 4, 1); //Create matrix for storing histograms

//Compile list of complete sets, and flag incomplete sets
currentRoot = "";
for(a=0; a<list.length; a++){
	//Make sure that the dataset is complete and follows standard naming convention
	root = replace(list[a], "-.*", "");
	if(matches(list[a], ".*-.*czi$") && !matches(root, currentRoot)){
		currentRoot = root;
		blankGel = "";
		blankSoln = "";
		gel = "";
		soln = "";
		for(b=0; b<list.length; b++){
			if(matches(list[b], root + "-blank-gel-fluorescence\\.czi$")) blankGel = list[b];
			if(matches(list[b], root + "-blank-soln-fluorescence\\.czi$")) blankSoln = list[b];
			if(matches(list[b], root + "-gel-fluorescence\\.czi$")) gel = list[b];
			if(matches(list[b], root + "-soln-fluorescence\\.czi$")) soln = list[b];
			if(!matches(blankGel, "") && !matches(blankSoln, "") && !matches(gel, "") && !matches(soln, "")){ //If dataset is complete, analyze set
				selectWindow("Normalized Ratio Histogram");
				while(nSlices < counter+1) run("Add Slice"); //Add slice to histogram stack if needed
				setSlice(nSlices);
				setMetadata("Label", root + ": " + histMin + " to " + histMax);
				quantifyFluorescence(blankGel, blankSoln, gel, soln, root); //This is where all of the analysis happens
				print("" + root + ": PASS");
				a=b;
				b=list.length;
				counter++;
			}
		}
		//If dataset is incomplete, report missing data in log file
		if(matches(blankGel, "")) print("Could not find: " + root + "-blank-gel-fluorescence.czi");
		if(matches(blankSoln, "")) print("Could not find: " + root + "-blank-soln-fluorescence.czi");
		if(matches(gel, "")) print("Could not find: " + root + "-gel-fluorescence.czi");
		if(matches(soln, "")) print("Could not find: " + root + "-soln-fluorescence.czi");
	}
}
setBatchMode("exit and display");
selectWindow("Results");
saveAs("Results", outDir + "Fluorescence quantification.csv");	
selectWindow("Normalized Ratio Histogram");
saveAs("Tiff", outDir + "Normalized Fluorescence Ratio Histograms.tiff");
selectWindow("Log");
saveAs("Text", outDir + "Processed file log.txt");


function quantifyFluorescence(blankGelStack, blankSolnStack, gelStack, solnStack, root){
	quantifyGelStack(gelStack);
	quantifySolnStack(solnStack);
	quantifyBlanks(blankGel, blankSoln);

	//Subtract respective blanks from gel or solution and report values
	imageCalculator("Subtract create 32-bit stack", "Soln", "SolnBlank");
	imageCalculator("Subtract create 32-bit stack", "Gel", "GelBlank");

	selectWindow("Result of Gel");
	Stack.getStatistics(voxelCount, mean, min, max, stdDev)
	row = nResults - 1;
	setResult("Gel minus blank", row, mean);

	selectWindow("Result of Soln");
	Stack.getStatistics(voxelCount, mean, min, max, stdDev)
	setResult("Soln minus blank", row, mean);
	
	// Calculate the partition coefficient, K
	imageCalculator("Divide create 32-bit stack", "Result of Gel","Result of Soln");
	
	selectWindow("Result of Result of Gel");
	Stack.getStatistics(voxelCount, Kmean, min, max, stdDev)
	setResult("K from soln above gel", row, Kmean);
	close("Result of Result of Gel");

	// Ali added 181006 to measure ratio of gel to solution from the coverslip
	imageCalculator("Subtract create 32-bit stack", "SolnCover", "SolnBlank");
	imageCalculator("Divide create 32-bit stack", "Result of Gel","Result of SolnCover");
	
	selectWindow("Result of Result of Gel");
	Stack.getStatistics(voxelCount, Kmean, min, max, stdDev)
	setResult("K from soln above coverslip", row, Kmean);

	// Close all open windows
	close("GelBlank");
	close("SolnBlank");
	close("Gel");	
	close("Soln");
	close("SolnCover");
	close("Result of Soln");
	close("Result of Gel");
	close("Result of Result of Gel");
	close("Result of SolnCover");
}

function quantifyBlanks(gel, soln){
	updateResults();
	for(a=0; a<2; a++){
		//Open the stack
		if(a){
			i=soln; 
			id = "Soln";
		}
		else{
			i=gel;
			id = "Gel";
		}
		row = nResults-1;
		Ext.openImagePlus(dir + i);
		selectWindow(i);
		slices = nSlices;

		offset = 36; // This corresponds to the slice that is 36 * 0.7 um = 25 um above starting point, since blank measurements are always started 20 um below the surface of the coverslip

	
		//Generate a Z-stack intensity profile
		selectWindow(i);
		analysisSlices = getResult("Gel Analysis Slices", row);
		start = offset;
		end = start + analysisSlices-1;
		run("Make Substack...", "  slices=" + start + "-" + end);
		selectWindow("Substack (" + start + "-" + end + ")");
		rename(id+"Blank");
		run("Median 3D...", "x=" + median + " y=" + median + " z=" + median); //Remove outlier pixels isotropicallyfileArray[a-offset] = mean;
		Stack.getStatistics(dummy, mean, dummy, dummy, dummy)
	
		setResult(id + "- blank", row, mean);

		//Calculate the ratio between the two stacks
		imageCalculator("Divide create 32-bit stack", id+"Blank","Soln");	
	
		//Get histogram of ratio
		selectWindow("Result of" + " " + id + "Blank");
		sumCounts = newArray(nBins);
		countSum = 0;
		for(c=1; c<=nSlices; c++){
			setSlice(c);
			getHistogram(values, counts, nBins, histMin, histMax);
			for(b=0; b<counts.length; b++){
				sumCounts[b] += counts[b];
				countSum += counts[b];
			}
		}
		//Export histogram 
		selectWindow("Normalized Ratio Histogram");
		for(b=0; b<nBins; b++){
			pixelValue = sumCounts[b]/countSum;
			if(isNaN(pixelValue)) pixelValue = 0;
			setPixel(b, 2+a, sumCounts[b]/countSum);
		}
	
		//Export results to results table
		selectWindow("Result of" + " " + id + "Blank");
		
		getVoxelSize(dummy, dummy, depth, unit);
		Stack.getStatistics(voxelCount, mean, min, max, stdDev)
		setResult("Blank  - " + id + " Ratio Mean", row, mean);
		setResult("Blank  - " + id + " Ratio Std Dev", row, stdDev);
		setResult("Blank  - " + id + " Ratio Min", row, min);
		setResult("Blank  - " + id + " Ratio Max", row, max);
		close("Result of" + " " + id + "Blank");
		close(i);
		//close(id+"Blank");
	}
	//close("Soln");
}


function quantifySolnStack(i){
	//Open the soln stack
	Ext.openImagePlus(dir + i);
	slices = nSlices;
	zProfileArray = newArray(slices);
	row = nResults-1;

	//Find the solution starting slice
	//Generate a Z-stack intensity profile
	selectWindow(i);
	for(a=0; a<slices; a++){
		setSlice(a+1);
		getStatistics(dummy, mean);
		zProfileArray[a] = mean;
	}

	// Ali added 181004
	//create duplicate image and the apply low pass filter to eliminate noise (gaussian blur)
	run("Duplicate...", "duplicate");
	run("Gaussian Blur 3D...", "x=1 y=1 z=1");

	zProfileArray_blur = newArray(slices);
	
	//Measure the Z profile of the gel stack
	for(a=0; a<slices; a++){
		setSlice(a+1);
		getStatistics(dummy, mean);
		zProfileArray_blur[a] = mean;
	}
	close();

	//take absolute derivative of intensity profile
	dzProfileArray_blur = newArray(slices);
	for(a=1; a<slices; a++){
		dzProfileArray_blur[a-1] = abs(zProfileArray_blur[a] - zProfileArray_blur[a-1]);
	}

	Array.getStatistics(dzProfileArray_blur, min, max, dummy, dummy);
	
	//Search for peak in derivative stack (glass->soln);
	dzPeakArray = newArray(0);
	dzPeakArray = Array.findMaxima(dzProfileArray_blur, max-min); 

	// Since there is no longer a glass slide on top, make the last slice the end of the solution region
	dzPeakArray = Array.concat(dzPeakArray,slices);
	
	solnArray = Array.slice(dzProfileArray_blur,dzPeakArray[0],dzPeakArray[1]);
	
	solnEndsArray = removeTransitions(solnArray, dzPeakArray[0]);
	//print(zProfileArray_blur[solnEndsArray[0]],zProfileArray_blur[solnEndsArray[1]]);

	gelEnd = getResult("Gel ending slice",row);

	solnStartCoverslip = solnEndsArray[0];

	//Make a substack of just the solution layer above the coverslip of equal thickness to prior analysis of the gel image, and compute the solution ratio
	selectWindow(i);

	endCoverslip = getResult("Gel Analysis Slices", row) + solnStartCoverslip - 1;
	setResult("Only Soln coverslip analysis  begin", row, solnStartCoverslip);
	setResult("Only Soln coverslip analysis end", row, endCoverslip);

	run("Make Substack...", "  slices=" + solnStartCoverslip + "-" + endCoverslip);
	selectWindow("Substack (" + solnStartCoverslip + "-" + endCoverslip + ")");
	rename("SolnCover");
	run("Median 3D...", "x=" + median + " y=" + median + " z=" + median); //Remove outlier pixels isotropically
	Stack.getStatistics(dummy, solnMeanCoverslip, dummy, dummy, dummy);
	setResult("Soln_coverlsip only value",row,solnMeanCoverslip);
	
	//Calculate the ratio between the two stacks (= image1/image2)
	imageCalculator("Divide create 32-bit stack", "Soln","SolnCover");	

	//Get histogram of ratio
	selectWindow("Result of Soln");
	sumCounts = newArray(nBins);
	countSum = 0;
	for(a=1; a<=nSlices; a++){
		setSlice(a);
		getHistogram(values, counts, nBins, histMin, histMax);
		for(b=0; b<counts.length; b++){
			sumCounts[b] += counts[b];
			countSum += counts[b];
		}
	}
	//Export histogram 
	selectWindow("Normalized Ratio Histogram");
	for(a=0; a<nBins; a++) setPixel(a, 1, sumCounts[a]/countSum);

	//Export results to results table
	selectWindow("Result of Soln");
	getVoxelSize(dummy, dummy, depth, unit);
	Stack.getStatistics(voxelCount, mean, min, max, stdDev)
	setResult("Soln above gel:Soln above coverslip Ratio Mean", row, mean);
	setResult("Soln above gel:Soln above coverslip Ratio Std Dev", row, stdDev);
	setResult("Soln above gel:Soln above coverslip Ratio Min", row, min);
	setResult("Soln above gel:Soln above coverslip Ratio Max", row, max);
	close("Result of Soln");
	//close(i);
	//close("Soln2");

	// Make a substack of just the solution layer at the same height as that used for the solution above the gel from the gel image
	solnStartGelMatch = solnEndsArray[0]+gelEnd;
	endSolnGel = getResult("Gel Analysis Slices", row) + solnStartGelMatch - 1;
	//print(solnStartGelMatch, endSolnGel);
	setResult("Only Soln same height as gel analysis  begin", row, solnStartGelMatch);
	setResult("Only Soln same height as gel end", row, endSolnGel);
	selectWindow(i);
	run("Make Substack...", "  slices=" + solnStartGelMatch + "-" + endSolnGel);
	selectWindow("Substack (" + solnStartGelMatch + "-" + endSolnGel + ")");
	rename("SolnControl");
	run("Median 3D...", "x=" + median + " y=" + median + " z=" + median); //Remove outlier pixels isotropically
	Stack.getStatistics(dummy, solnMeanGel, dummy, dummy, dummy);
	setResult("Soln only value at gel height",row,solnMeanGel);

	//Calculate the ratio between the two stacks (= image1/image2)
	imageCalculator("Divide create 32-bit stack", "SolnControl","SolnCover");
	selectWindow("Result of SolnControl");
	getVoxelSize(dummy, dummy, depth, unit);
	Stack.getStatistics(voxelCount, mean, min, max, stdDev)
	setResult("Soln control:Soln coverslip Ratio Mean", row, mean);
	setResult("Soln control:Soln coverslip Ratio Std Dev", row, stdDev);
	setResult("Soln control:Soln coverslip Ratio Min", row, min);
	setResult("Soln control:Soln coverslip Ratio Max", row, max);
	close("SolnControl");
	close("Result of SolnControl");
	close(i);
}
function quantifyGelStack(i){
	//Open the gel stack
	Ext.openImagePlus(dir + i);
	slices = nSlices;
	zProfileArray_orig = newArray(slices);
	zProfileArray_blur = newArray(slices);
	row = nResults;
	setResult("Sample", row, root);

	//Measure the Z profile of the gel stack and save to "original" array
	for(a=0; a<slices; a++){
		setSlice(a+1);
		getStatistics(dummy, mean);
		zProfileArray_orig[a] = mean;
	}

	//AS 180923: create duplicate image and the apply low pass filter to eliminate noise (gaussian blur) for finding the region ends
	run("Duplicate...", "duplicate");
	run("Gaussian Blur 3D...", "x=1 y=1 z=1");
	
	//Measure the Z profile of the gel stack and save to "blurred" array
	for(a=0; a<slices; a++){
		setSlice(a+1);
		getStatistics(dummy, mean);
		zProfileArray_blur[a] = mean;
	}

	//Close the Gaussian blurred image
	close();

	//take absolute derivative of intensity profile
	dzProfileArray_orig = newArray(slices);
	dzProfileArray_blur = newArray(slices);
	for(a=1; a<slices; a++){
		dzProfileArray_orig[a-1] = abs(zProfileArray_orig[a] - zProfileArray_orig[a-1]);
		dzProfileArray_blur[a-1] = abs(zProfileArray_blur[a] - zProfileArray_blur[a-1]);
	}
	//Array.show(dzProfileArray_blur);
	//print(dzProfileArray_blur[0]);
	
	//Search for 2 peaks in blurred derivative stack (glass->gel, gel->soln);
	dzPeakArray = newArray(0);
	Array.getStatistics(dzProfileArray_blur, min, max, dummy, dummy);

	// max-min is the tolerance (minimum amplitude difference needed to separate two peaks)
	// slowly decrease the tolerance until find the 2 peaks
	// findMaxima returns the peak positions in the array (not values). dzPeakArray contains the location in the profile array of the derivative peaks, sorted in desceinding strength (written into findMaxima function)
	while(dzPeakArray.length < 2){
		dzPeakArray = Array.findMaxima(dzProfileArray_blur, max-min); 
		max--;
	}

	// Since there is no longer a glass slide on top, make the last slice the end of the solution region
	dzPeakArray = Array.concat(dzPeakArray,slices);

	// We want the peak array sorted from first to last location
	Array.sort(dzPeakArray);
	//Array.show(dzPeakArray);

	// Generate gel and soln derivative arrays from the gel image by extracting the regions between the peaks 
	gelArray = Array.slice(dzProfileArray_blur,dzPeakArray[0],dzPeakArray[1]); //array.slice extracts part of an array and returns it 
	solnArray = Array.slice(dzProfileArray_blur,dzPeakArray[1],dzPeakArray[2]);

	//Array.show(gelArray);

	// Currently in the transition region of the profiles - need to find the correct edges for quantification
	gelEndsArray = removeTransitions(gelArray, dzPeakArray[0]);
	solnEndsArray = removeTransitions(solnArray, dzPeakArray[1]);

	//print(zProfileArray_blur[gelEndsArray[0]], zProfileArray_blur[gelEndsArray[1]], zProfileArray_blur[solnEndsArray[0]],zProfileArray_blur[solnEndsArray[1]]);
	
	//Create a gel and a solution substack of equal thickness (to maintain equal information content): depending on whether gel or solution substack is larger to start with
	if((gelEndsArray[1] - gelEndsArray[0]) <= (solnEndsArray[1] - solnEndsArray[0])){
		selectWindow(i);
		run("Make Substack...", "  slices=" + gelEndsArray[0] + "-" + gelEndsArray[1]); // Make gel substack from identified endpoints 
		selectWindow("Substack (" + gelEndsArray[0] + "-" + gelEndsArray[1] + ")");
		rename("Gel");
		run("Median 3D...", "x=" + median + " y=" + median + " z=" + median); //Remove outlier pixels isotropically
		Stack.getStatistics(dummy, ingelMean, dummy, dummy, dummy);
		nGelSlices = nSlices;
		selectWindow(i); // Make solution substack from identified endpoints, cutting off the end to the same length as the gel substack
		run("Make Substack...", "  slices=" + solnEndsArray[0] + "-" + (solnEndsArray[0] + nGelSlices - 1));
		selectWindow("Substack (" + solnEndsArray[0] + "-" + (solnEndsArray[0] + nGelSlices - 1) + ")");
		rename("Soln");
		run("Median 3D...", "x=" + median + " y=" + median + " z=" + median); //Remove outlier pixels isotropically
		Stack.getStatistics(dummy, insolnMean, dummy, dummy, dummy);
		setResult("Limiting thickness", row, "Gel"); // Report that the gel was the limited thickness
		setResult("Gel starting slice",row,gelEndsArray[0]);
		setResult("Gel ending slice",row,gelEndsArray[1]);
		setResult("Soln starting slice",row,solnEndsArray[0]);
		setResult("Soln ending slice",row,(solnEndsArray[0] + nGelSlices - 1));
	}
	else{
		selectWindow(i);
		run("Make Substack...", "  slices=" + solnEndsArray[0] + "-" + solnEndsArray[1]);
		selectWindow("Substack (" + solnEndsArray[0] + "-" + solnEndsArray[1] + ")");
		rename("Soln");
		run("Median 3D...", "x=" + median + " y=" + median + " z=" + median); //Remove outlier pixels isotropically
		Stack.getStatistics(dummy, insolnMean, dummy, dummy, dummy);
		nSolnSlices = nSlices;
		selectWindow(i);
		run("Make Substack...", "  slices=" + gelEndsArray[0] + "-" + (gelEndsArray[0] + nSolnSlices - 1));
		selectWindow("Substack (" + gelEndsArray[0] + "-" + (gelEndsArray[0] + nSolnSlices - 1) + ")");
		rename("Gel");
		run("Median 3D...", "x=" + median + " y=" + median + " z=" + median); //Remove outlier pixels isotropically
		Stack.getStatistics(dummy, ingelMean, dummy, dummy, dummy);
		setResult("Limiting thickness", row, "Solution"); 
		setResult("Gel starting slice",row,gelEndsArray[0]);
		setResult("Gel ending slice",row,(gelEndsArray[0] + nSolnSlices - 1));
		setResult("Soln starting slice",row,solnEndsArray[0]);
		setResult("Soln ending slice",row,solnEndsArray[1]);
	}

	setResult("Gel value",row,ingelMean);
	setResult("Soln value",row,insolnMean);
	
	//Calculate the ratio between the two stacks
	imageCalculator("Divide create 32-bit stack", "Gel","Soln");	

	//Get histogram of ratio
	selectWindow("Result of Gel");
	sumCounts = newArray(nBins);
	countSum = 0;
	for(a=1; a<=nSlices; a++){
		setSlice(a);
		getHistogram(values, counts, nBins, histMin, histMax);
		for(b=0; b<counts.length; b++){
			sumCounts[b] += counts[b];
			countSum += counts[b];
		}
	}
	//Export histogram 
	selectWindow("Normalized Ratio Histogram");
	for(a=0; a<nBins; a++) setPixel(a, 0, sumCounts[a]/countSum);

	//Export results to results table
	selectWindow("Result of Gel");
	getVoxelSize(dummy, dummy, depth, unit);
	Stack.getStatistics(voxelCount, mean, min, max, stdDev)
	setResult("Gel Thickness (" + unit + ")", row, (dzPeakArray[1]-dzPeakArray[0])*depth);
	setResult("Gel Analysis Thickness (" + unit + ")", row, nSlices*depth);
	setResult("Gel Analysis Slices", row, nSlices);
	setResult("Gel Ratio Mean", row, mean);
	setResult("Gel Ratio Std Dev", row, stdDev);
	setResult("Gel Ratio Min", row, min);
	setResult("Gel Ratio Max", row, max);
	setResult("Gel Ratio Counts (nVoxels)", row, voxelCount);
	close("Result of Gel");
	close(i);
	//close("Gel");
}

//Remove transition regions from arrays
function removeTransitions(array, offset){
	start = 0;
	end = array.length;

// 180923: Ali added looking for second derivative. Create array to store second derivative of the correct length
	dzdzProfileArray = newArray(array.length);

	// Calculate second derivative
	for(a=1; a<array.length; a++){
		dzdzProfileArray[a-1] = array[a] - array[a-1];
	}

	// Find end of transition region from the very left by where the second derivative changes sign. Since we know we are going from glass into gel or gel into solution, the second derivative should go from negative to positive
	for(b=0; b<dzdzProfileArray.length; b++){
		if(dzdzProfileArray[b] >0){
			start = b;
			b=dzdzProfileArray.length;	
		}
	}

	// Find end of transition region from the right by where the second derivative changes sign. Now, since going from solution to gel, the second derivative should go from positive to negative
	// If the derivative doesn't change signs, leave it as the last frame (true for solution)
	for(b=dzdzProfileArray.length-1; b>start; b--){
		if(dzdzProfileArray[b] <0){
			end = b;
			b=start;	
		}
	}
	//print(dzdzProfileArray[start],dzdzProfileArray[end]);

	// Return the derivative array starting at the correct index from the beginning of the transition region identified previously
	returnArray2 = newArray(offset + start + 1, offset + end - 1);
	return returnArray2;
}


