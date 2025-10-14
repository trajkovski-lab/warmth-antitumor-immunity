//This macro will let us check each image one by one and correct manually the segmentations made automatically by CellPose if necessary

extension=".png";


dir = getDirectory("Select a directory containing one or several lif files.");
files = getFileList(dir);

for(f=0; f<files.length; f++) {
	if(endsWith(files[f], extension)) {
		open(dir + files[f]);
		image=getTitle();
		
		if(!File.exists(dir + image + "_RESULTSFINALCORRECTED.txt")){
			run("RGB Color");
			run("8-bit");
			roiManager("Open", dir + image + "_ch1_allcrops0.zip"); //adjust here depending on channel segmented
			roiManager("Show All");
			waitForUser("Check");
			roiManager("Deselect");
			roiManager("Measure");
			saveAs("Results", dir + image + "_RESULTSFINALCORRECTED.txt");
			roiManager("Deselect");
			roiManager("Save", dir + image + "_corrected.zip");
			selectWindow("Results"); run("Close");
			selectWindow("ROI Manager"); run("Close");
			run("Close All");
		} else {
			run("Close All");
		}
		
	}
}
