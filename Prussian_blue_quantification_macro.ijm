// Select a folder of images
MyDir = getDirectory("mypath");
inputDir = MyDir + "Images/";

// Create output directories
resultDir = MyDir + "Results/";
maskDir = MyDir + "Masks/";
if (!File.exists(resultDir)) File.makeDirectory(resultDir);
if (!File.exists(maskDir)) File.makeDirectory(maskDir);


// Process all files
list = getFileList(inputDir);
for (i = 0; i < list.length; i++) {
    filename = list[i];
    fullPath = inputDir + filename;

    // Open image
   run("Bio-Formats Importer", "open=[" + fullPath + "] autoscale color_mode=Composite open_all_series series=3 windowless=true");

    // Image processing pipeline
    run("RGB Color");
    run("Colour Deconvolution2", "vectors=[FastRed FastBlue DAB] output=RGB_Keep_C2 simulated cross hide");
    run("8-bit");
    run("Auto Threshold", "method=IsoData white");
    setOption("BlackBackground", true);
    run("Make Binary");
    run("Convert to Mask");
    run("Erode");
    run("Watershed");
    run("Find Maxima...", "prominence=100 strict output=[Point Selection]");
    
    // Particle analysis with proper row handling
    run("Analyze Particles...", "size=50-Infinity pixel display clear summarize overlay");
    
    // Get current row count safely
    rowCount = nResults();
    
    // Save individual results if particles found
    if (isOpen("Results")) {
        selectWindow("Results");
        imgName = getTitle();
        imgBase = substring(imgName, 0, lastIndexOf(filename, "."));
        saveAs("Results", resultDir+imgBase+"_Results.csv");
    }
    
     // Add zero entry if no particles found (using proper row index)
    if (rowCount == 0) {
        setResult("Image", rowCount, filename);
        setResult("Count", rowCount, 0);
        setResult("Area", rowCount, 0);
        updateResults();
    }

    // Save mask
    saveAs("Tiff", maskDir+imgBase+"_Mask.tif");
    while (nImages > 0)
    close();
}

// Save final summary
selectWindow("Summary");
saveAs("Results", resultDir+"Summary.csv");

run("Close All");
