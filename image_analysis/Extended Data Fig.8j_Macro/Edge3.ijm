//This macro lets you manually draw the tumor edge, and the proximal and distal regions, to analyse presegmented lipid droplets separately


dir=getDirectory("image");
title=getTitle();

roiManager("deselect");
roiManager("save", dir+title+"_LDs.zip"); 
nld=roiManager("count");


roiManager("deselect");


waitForUser("Draw Tumour Edge");
roiManager("add");
edge=roiManager("count")-1;

roiManager("Select", edge);
run("Interpolate", "interval=1 smooth");
roiManager("add");
roiManager("Select", edge);
roiManager("save", dir+ title + "_EdgeDrawn.roi");
roiManager("Delete");
roiManager("Select", edge);
roiManager("save", dir+ title + "_Edge.roi");

roiManager("Select", edge);
Roi.getCoordinates(xpoints, ypoints);

//Draw two categories and use them separately

waitForUser("Draw two regions (tumour first) and press T key after each of them");

nld1=nld+1;
nld2=nld+2;

roiManager("deselect");






aa=nld;
bb=nld+1;
cc=nld+2;
dd=nld+3;

allplus0=Array.getSequence(aa); // fom 0 to nld-1. a sequence of size nld
allplus1=Array.concat(allplus0, bb); //we add here the first group selection
allplus2=Array.concat(allplus0, cc); //we add here the second group selection
roiManager("Select", allplus1);

roiManager("Combine");
wait(3000);
roiManager("Add"); //this is nld+3.

roiManager("Select", allplus2);
roiManager("Combine");
wait(3000);
roiManager("Add"); //this is nld+4. 

todelete=Array.getSequence(dd);
roiManager("Deselect");
roiManager("Select", todelete);
wait(2000);
roiManager("Delete");

roiManager("Deselect");
roiManager("Save", dir+title+"_RegionComposites.zip");


//Now we have only two ROI

roiManager("Select", 0);
roiManager("Split");
lastroi2=roiManager("Count"); 


roiManager("Select", 1);
roiManager("Split");
lastroi1=roiManager("Count"); 
//group 1 goes from ROI lastroi2+1 to lastroi1
toendof1=Array.getSequence(lastroi2);
adipos1=Array.deleteIndex(toendof1, 0);
adipos1=Array.deleteIndex(adipos1, 0);


roiManager("Deselect");

roiManager("Save", dir+title+"_regionsAndLDs.zip");
wait(3000);
roiManager("Deselect");

//we should have now a ROI Manager with the two composites and all the adipocytes again, but now ordered by group




//now lets loop through the adipocytes once for each group
total=roiManager("Count")-2

distances=newArray(total);
nearxcoords=newArray(total);
nearycoords=newArray(total);
id=newArray(total);
xLD=newArray(total);
yLD=newArray(total);
class=newArray(total);
area=newArray(total);
circularity=newArray(total);



roiManager("Deselect");
roiManager("Measure"); 


//GROUP 2
for(i=2; i<lastroi2; i++) {
	roiManager("Select", i);
	X=getResult("X", i);
	Y=getResult("Y", i);
	LDarea=getResult("Area", i);
	circ=getResult("Circ.", i);
	
	if (LDarea < 20000) {
		mindist=100000000000000;
		for(j=0; j<xpoints.length; j++) {
			eX=xpoints[j];
			eY=ypoints[j];
			distance=sqrt(pow((eX-X),2)+pow((eY-Y),2));
			if (distance < mindist) {
				mindist=distance;
				nearX=eX;
				nearY=eY;
			}	
		}
		distances[i-2]=mindist;
		nearxcoords[i-2]=nearX;
		nearycoords[i-2]=nearY;
		id[i-2]=i+1;
		xLD[i-2]=X;
		yLD[i-2]=Y;
		class[i-2]="OUTSIDE";
		area[i-2]=LDarea;
		circularity[i-2]=circ;
	} else {
		distances[i-2]=NaN;
		nearxcoords[i-2]=NaN;
		nearycoords[i-2]=NaN;
		id[i-2]=NaN;
		xLD[i-2]=NaN;
		yLD[i-2]=NaN;
		class[i-2]=NaN;
		area[i-2]=NaN;
		circularity[i-2]=NaN;
	}
} 



//GROUP 1
for(i=lastroi2+1; i<lastroi1; i++) {
	roiManager("Select", i);
	X=getResult("X", i);
	Y=getResult("Y", i);
	LDarea=getResult("Area", i);
	circ=getResult("Circ.", i);
	
	if (LDarea < 20000) { //when some adipocytes are excluded from both classess, the region may be included in the ROI list
		mindist=100000000000000;
		for(j=0; j<xpoints.length; j++) {
			eX=xpoints[j];
			eY=ypoints[j];
			distance=sqrt(pow((eX-X),2)+pow((eY-Y),2));
			if (distance < mindist) {
				mindist=distance;
				nearX=eX;
				nearY=eY;
			}	
		}
		distances[i-3]=mindist;
		nearxcoords[i-3]=nearX;
		nearycoords[i-3]=nearY;
		id[i-3]=i+1;
		xLD[i-3]=X;
		yLD[i-3]=Y;
		class[i-3]="INSIDE";
		area[i-3]=LDarea;
		circularity[i-3]=circ;
	} else {
		distances[i-3]=NaN;
		nearxcoords[i-3]=NaN;
		nearycoords[i-3]=NaN;
		id[i-3]=NaN;
		xLD[i-3]=NaN;
		yLD[i-3]=NaN;
		class[i-3]=NaN;
		area[i-3]=NaN;
		circularity[i-3]=NaN;
	}
}


//////////////////////////////////////
///Check if any adipocyte repeated
//That would mean it was not included in any class

location=newArray(total);
for (i=0; i<distances.length; i++){
	ldx1=xLD[i];
	ldy1=yLD[i];
	classi=class[i];
	
	counter=0;
	for (j=0; j<distances.length; j++){
		ldx2=xLD[j];
		ldy2=yLD[j];
		if (ldx1==ldx2 && ldy1==ldy2){
			counter=counter+1;
		}
	}
	
	if (counter >1) {
		location[i] = "EXCLUDED";
	} else {
		location[i] =classi;
	}
}









///////////////////////////////////////






Array.show(id, distances, nearxcoords, nearycoords, xLD, yLD, class,circularity, area, location);
selectWindow("Arrays");





saveAs("results", dir + title + "_Distances.txt");



selectWindow(title + "_Distances.txt");
run("Close");
selectWindow("Results");
saveAs("results", dir + title + "_ResultsTable.txt");
run("Close");

roiManager("Deselect");
roiManager("Select", toendof1);

roiManager("Delete");
roiManager("Select", 0); 
roiManager("Delete");
last=roiManager("Count")-1;


roiManager("save", dir+title+"_insideLDs.zip");
selectWindow("ROI Manager");
run("Close");


roiManager("open", dir+title+"_RegionsAndLDs.zip");
all=roiManager("count");
select0=Array.getSequence(all);
select1=Array.slice(select0, lastroi2, all);

roiManager("Deselect");
roiManager("Select", select1); 

roiManager("Delete");
roiManager("Deselect");
select2=Array.getSequence(2);
roiManager("Select", select2);
roiManager("Delete");

roiManager("Deselect");
roiManager("save", dir+title+"_outsideLDs.zip");

selectWindow("ROI Manager");
run("Close");




selectWindow(title);
run("Close");









