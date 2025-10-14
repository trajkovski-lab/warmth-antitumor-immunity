//This macro corrects segmentations produced by MRI Adipocyte Tools based on circularity and intensity

run("Set Measurements...", "area mean min centroid shape redirect=None decimal=3");
nld=roiManager("count");
print(nld);

roiManager("deselect");
roiManager("measure");


n=0
for (i=0; i<nld; i++) {
	roiManager("select", n);
	
	if((round(n/500) - n/500) == 0) {
		print(n);
	}
	
	selectWindow("Results");
	//wait(1);
	intensity=getResult("Mean", i);
	circularity=getResult("Circ.", i);

	
	if (intensity < 200 || circularity < 0.15) {
		roiManager("select", n);
		roiManager("delete");
		n=n-1;
	}
	n=n+1;	
}

nld=roiManager("count");
print(nld);


