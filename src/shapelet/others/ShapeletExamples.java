package shapelet.others;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;
import com.lw.shapelet.*;
import com.opencsv.CSVWriter;
import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.evaluation.NominalPrediction;
import weka.classifiers.trees.RandomForest;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.SerializationHelper;
import weka.core.Utils;
import weka.filters.timeseries.shapelet_transforms.*;

public class ShapeletExamples {

	public static ShapeletTransform st;
	public static ArrayList<Shapelet> clusteredShapelets;
	public static ArrayList<Shapelet> max_clusteredShapelets;

	public static Instances basicTransformExample(Instances train) {
		st = new ShapeletTransform();
		System.out.println("classsssssssssssss = "+train.instance(0).classValue());
		// Let m=train.numAttributes()-1 (series length)
		// Let n= train.numInstances() (number of series)
		int nosShapelets = (train.numAttributes() - 1) * train.numInstances() / 5;
		if (nosShapelets < ShapeletTransform.DEFAULT_NUMSHAPELETS)
			nosShapelets = ShapeletTransform.DEFAULT_NUMSHAPELETS;
		st.setNumberOfShapelets(nosShapelets);
		int minLength = 5;
		int maxLength = (train.numAttributes() - 1) / 2;
		System.out.println("maxlength = " + maxLength);
		if (maxLength < ShapeletTransform.DEFAULT_MINSHAPELETLENGTH)
			maxLength = ShapeletTransform.DEFAULT_MINSHAPELETLENGTH;
		st.setShapeletMinAndMax(minLength, maxLength);
		st.setQualityMeasure(QualityMeasures.ShapeletQualityChoice.F_STAT);
		st.setLogOutputFile("ShapeletExampleLog.csv");
		Instances shapeletT = null;
		try {
			shapeletT = st.process(train);
		} catch (Exception ex) {
			System.out.println("Error performing the shapelet transform" + ex);
			ex.printStackTrace();
			System.exit(0);
		}
		return shapeletT;
	}

	public static Instances clusteredShapeletTransformExample(Instances train) {
		Instances shapeletT = null;

		int nosShapelets = (train.numAttributes() - 1) * train.numInstances() / 50;
		// System.out.println("noshapelets = "+nosShapelets);
		ClusteredShapeletTransform cst = new ClusteredShapeletTransform(st, nosShapelets);
		System.out.println(" Clustering down to " + nosShapelets + " Shapelets");
		System.out.println(" From " + st.getNumberOfShapelets() + " Shapelets");
		// ShapeletTransform st = new ShapeletTransform();
		try {
			shapeletT = cst.process(train);
			// Instances output = cst.calDistance(train,shapelets);
			// return output;
		} catch (Exception ex) {
			System.out.println("Error performing the shapelet clustering" + ex);

			ex.printStackTrace();
			System.exit(0);
		}
		clusteredShapelets = cst.clusteredShapelets;
		return shapeletT;

	}

	public static void initializeShapelet(ShapeletTransform s, Instances train) {
		s.setNumberOfShapelets(1);
		int minLength = 15;
		int maxLength = 36;
		s.setShapeletMinAndMax(minLength, maxLength);
		s.setQualityMeasure(QualityMeasures.ShapeletQualityChoice.F_STAT);
		s.supressOutput();
		s.turnOffLog();
	}
	public static void main(String[] args) throws Exception {
		Instances alldata = null, train = null, test = null;
		FileReader r;
		for (int a = 0; a < 6; a++) {
			try {
				r = new FileReader("wrist_" + a + "_column.csv");
				alldata = new Instances(r);
				alldata.setClassIndex(alldata.numAttributes() - 1);
			}
			
			catch (Exception e) {
				System.out.println("Unable to load data. Exception thrown =" + e);
				System.exit(0);
			}
			System.out.println("------------- = "+alldata.instance(0).classValue());
			//System.out.println(alldata.instance(0).classValue());
            System.out.println(alldata.classAttribute().value((int) alldata.instance(0).classValue()));
			System.out.println(alldata.instance(0).toString().split(",")[alldata.instance(0).toString().split(",").length-1]);
			Instances shapeletT;
			Instance instancd;
			Instances shapeletCluster;
			int seed = 0;
			int folds = 3;
			int runs = 3;
			int sum_j48 = 0;
			int sum_randomforest = 0;
			int sum_adtree = 0;
			double max = 0;
		//	for (int n = 0; n < 2; n++) {
				seed += 1;
				Random rand = new Random(seed);
				Instances randData = new Instances(alldata);
				randData.randomize(rand);
				if (randData.classAttribute().isNominal())
					randData.stratify(folds);

				//for (int i = 0; i < folds; i++) {				
					train = randData.trainCV(folds, 1);
					test = randData.testCV(folds, 1);
                    System.out.println(train.numInstances());
                    System.out.println(test.numInstances());
                    shapeletT = basicTransformExample(train);
					shapeletCluster = clusteredShapeletTransformExample(train);
					writeToFile("classin.csv", shapeletCluster);
					int noClust = shapeletCluster.numAttributes() - 1;
					Instances instancesTest = calDistance(noClust, test, clusteredShapelets);
					System.out.println("=======================\n instances test = \n");
					writeToFile("classtest.csv", instancesTest);

					// * 璁粌妯″瀷锛岄娴�
					Instances train1 = null, test1 = null;
					FileReader r1;
					try {
						r1 = new FileReader("classin.csv");
						train1 = new Instances(r1);
						train1.setClassIndex(train1.numAttributes() - 1);
						// System.out.println("read train1 end===========");
						r1 = new FileReader("classtest.csv");
						test1 = new Instances(r1);
						test1.setClassIndex(test1.numAttributes() - 1);
						// System.out.println("read test1 end=========");
					} catch (Exception e) {
						System.out.println("Unable to load data. Exception thrown =" + e);
						System.exit(0);
					}
					// FastVector predictions_j48 = new FastVector();
					// Evaluation evaluation_j48 = new Evaluation(train1);

					FastVector predictions_RandomForest = new FastVector();
					Evaluation evaluation_RandomForest = new Evaluation(train1);

					// FastVector predictions_adtree = new FastVector();
					// Evaluation evaluation_adtree = new Evaluation(train1);

					// Classifier model_j48 = new J48();
					Classifier model_RandomForest = new RandomForest();
					// Classifier model_adtree = new ADTree();

					// model_j48.buildClassifier(train1);
					model_RandomForest.buildClassifier(train1);
					SerializationHelper.write("LibSVM.model", model_RandomForest);
					// model_adtree.buildClassifier(train1);

					// evaluation_j48.evaluateModel(model_j48, test1);
					evaluation_RandomForest.evaluateModel(model_RandomForest, test1);

					// evaluation_adtree.evaluateModel(model_adtree, test1);

					// predictions_j48.appendElements(evaluation_j48.predictions());
					predictions_RandomForest.appendElements(evaluation_RandomForest.predictions());
					// predictions_adtree.appendElements(evaluation_adtree.predictions());

					// double accuracy_j48 = calculateAccuracy(predictions_j48);
					double accuracy_RandomForest = calculateAccuracy(predictions_RandomForest);
					// double accuracy_adtree =
					// calculateAccuracy(predictions_adtree);

					/*
					 * System.out.println("Accuracy of " + new
					 * J48().getClass().getSimpleName() + ": " +
					 * String.format("%.2f%%", accuracy_j48) +
					 * "\n---------------------------------");
					 */
					System.out.println("Accuracy of " + new RandomForest().getClass().getSimpleName() + ": "
							+ String.format("%.2f%%", accuracy_RandomForest) + "\n---------------------------------");
					// System.out.println("Accuracy of " + new
					// ADTree().getClass().getSimpleName() + ": "
					// + String.format("%.2f%%", accuracy_adtree)
					// + "\n---------------------------------");
					/*if (accuracy_RandomForest > max) {
						max = accuracy_RandomForest;
						max_clusteredShapelets = clusteredShapelets;
					}*/
					// sum_j48 += accuracy_j48;
				//	sum_randomforest += accuracy_RandomForest;
					// sum_randomforest += accuracy_RandomForest;
					// sum_adtree +=accuracy_adtree;

			//	}
		//	}
			System.out.println("random forest average accuracy = " + (accuracy_RandomForest));
			System.out.println("涓嬩竴涓�------------------------------------------------------");
			writeShapelet(clusteredShapelets, a);
		}

	}

	private static void writeShapelet(ArrayList<Shapelet> clusteredShapelets, int a) throws IOException {
		// TODO Auto-generated method stub
		// FileWriter writer = new FileWriter("shapelet.csv",true);
		CSVWriter writer = new CSVWriter(new FileWriter("shapelet.csv", true));
		for (int i = 0; i < clusteredShapelets.size(); i++) {
			double[] shapelet = clusteredShapelets.get(i).getContent();
			String[] shapelet_string = new String[shapelet.length + 2];
			// writer.writeNext(shapelet);
			for (int j = 0; j < shapelet.length; j++) {
				shapelet_string[j] = shapelet[j] + "";
				if (shapelet_string[j].equals("NAN")) {
					System.out.println("绗�" + i + "涓猻hapelet鐨勭" + j + "涓�兼槸绌�");
				}
			}
			shapelet_string[shapelet.length] = a + "";
			shapelet_string[shapelet.length + 1] = clusteredShapelets.get(i).getShapeletClass() + "";
			writer.writeNext(shapelet_string);

		}
		writer.close();
	}

	public static void writeToFile(String fileName, Instances shapeletT) throws IOException {
		StringBuffer text = new StringBuffer();
		text.append(shapeletT.ARFF_RELATION).append(" ").append(Utils.quote(shapeletT.relationName())).append("\n\n");
		for (int i = 0; i < shapeletT.numAttributes(); i++) {
			text.append(shapeletT.attribute(i)).append("\n");
		}
		text.append("\n").append(shapeletT.ARFF_DATA).append("\n");
		StringBuffer text1 = new StringBuffer();
		// System.out.println("shapeletT.numInstances() =
		// "+shapeletT.numInstances());
		for (int j = 0; j < shapeletT.numInstances(); j++) {
			text1.append(shapeletT.instance(j));
			if (j < shapeletT.numInstances() - 1) {
				text1.append('\n');
			}
		}

		text.append(text1.toString());
		FileWriter writer = new FileWriter(fileName);
		writer.write(text.toString());
		writer.close();
	}

	public static double calculateAccuracy(FastVector predictions) {
		double correct = 0;

		for (int i = 0; i < predictions.size(); i++) {
			NominalPrediction np = (NominalPrediction) predictions.elementAt(i);
			if (np.predicted() == np.actual()) {
				correct++;
			}
		}
		return 100 * correct / predictions.size();
	}

	public static Instances calDistance(int noClust, Instances data, ArrayList<Shapelet> shapeletT) throws Exception {

		Instances output = determineOutputFormat(noClust, data);
		for (int i = 0; i < data.numInstances(); i++) { // for each data
			Instance toAdd = new Instance(shapeletT.size() + 1);
			int shapeletNum = 0;
			for (Shapelet s : shapeletT) {
				double dist = ShapeletTransform.subsequenceDistance(s.content, data.instance(i));
				toAdd.setValue(shapeletNum++, dist);
			}
			toAdd.setValue(shapeletT.size(), data.instance(i).classValue());
			output.add(toAdd);
		}
		return output;
	}

	public static Instances determineOutputFormat(int noClust, Instances inputFormat) throws Exception {

		FastVector atts = new FastVector();
		String name;
		for (int i = 0; i < noClust; i++) {
			name = "CShapelet_" + i;
			atts.addElement(new Attribute(name));
		}
		Attribute target = inputFormat.attribute(inputFormat.classIndex());
		FastVector vals = new FastVector(target.numValues());
		for (int i = 0; i < target.numValues(); i++) {
			vals.addElement(target.value(i));
		}
		atts.addElement(new Attribute(inputFormat.attribute(inputFormat.classIndex()).name(), vals));
		Instances result = new Instances("CShapelets" + inputFormat.relationName(), atts, inputFormat.numInstances());
		result.setClassIndex(result.numAttributes() - 1);
		return result;
	}
	
	 /* public static void distanceOptimizations(Instances train){ 
		  Instances shapeletT=null; ShapeletTransform s1=new ShapeletTransform();
	  initializeShapelet(s1,train); // ShapeletTransformOnlineNorm s2=new
	  ShapeletTransformOnlineNorm(); // initializeShapelet(s2,train); //
	  ShapeletTransformDistCaching s3=new ShapeletTransformDistCaching(); //
	  initializeShapelet(s3,train); 
	  DecimalFormat df =new DecimalFormat("###.####"); long t1=0; long t2=0; double
	  time1,time2,time3; try { t1=System.currentTimeMillis();
	  shapeletT=s1.process(train); t2=System.currentTimeMillis();
	  time1=((t2-t1)/1000.0); t1=System.currentTimeMillis(); //
	  shapeletT=s2.process(train);
	  t2=System.currentTimeMillis();
	  time2=((t2-t1)/1000.0); 
	  t1=System.currentTimeMillis(); //
	  shapeletT=s3.process(train); 
	  t2=System.currentTimeMillis();
	  time3=((t2-t1)/1000.0); 
	  System.out.println("TIME (seconds)"); 
	  System.out.println("No Optimization\t Online Norm/Early Abandon\t Distance caching"
	  );
	  System.out.println(df.format(time1)+"\t\t\t"+df.format(time2)+"\t\t\t"+df.format(time3)); System.out.
	  println("TIME REDUCTION\t Online Norm/Early Abandon\t Distance caching");
	  System.out.println("\t\t\t"+(int)(100.0*time2/time1)+"% \t\t\t"+(int)(100.0*time3/time1)+"%"); System.out.
	  println("SPEED UP\t Online Norm/Early Abandon\t Distance caching");
	  System.out.println("\t\t\t"+df.format(time1/time2)+"\t\t\t"+df.format(
	  time1/time3)); } catch (Exception ex) {
	  System.out.println("Error performing the shapelet transform"+ex);
	  ex.printStackTrace(); System.exit(0); } } public static void
	  shapeletEarlyAbandons(Instances train){ 
		  //Time the speed up from early abandon of the four distance measures.
	  
	  //IG: 
		  ShapeletTransform[] s=new ShapeletTransform[4]; 
		  ShapeletTransform[] pruned=new ShapeletTransform[4]; for(int i=0;i<s.length;i++){ 
			  s[i]=new ShapeletTransformDistCaching(); 
			  pruned[i]=new ShapeletTransformDistCaching(); 
			  } 
		  for(ShapeletTransform s1:s){
	  initializeShapelet(s1,train); 
	  s1.setCandidatePruning(false); 
	  }
	  for(ShapeletTransform s1:pruned){ 
		  initializeShapelet(s1,train);
	  s1.setCandidatePruning(true); 
	  } 
	  QualityMeasures.ShapeletQualityChoice[] choices=QualityMeasures.ShapeletQualityChoice.values(); 
	  for(int i=0;i<s.length;i++){
		  s[i].setQualityMeasure(choices[i]);
	  pruned[i].setQualityMeasure(choices[i]);
	  } 
	  long t1,t2;
	  double time1,time2; 
	  DecimalFormat df =new DecimalFormat("###.####"); 
	  try {
	  for(int i=0;i<s.length;i++){ 
		  t1=System.currentTimeMillis();
	 s[i].process(train); 
	 t2=System.currentTimeMillis();
	  time1=((t2-t1)/1000.0); 
	  t1=System.currentTimeMillis();
	  pruned[i].process(train); 
	  t2=System.currentTimeMillis();
	  time2=((t2-t1)/1000.0);
	  System.out.println(" ********* QUALITY MEASURE ="+s[i].getQualityMeasure()+"  **********"); System.out.
	  println(" NO ABANDON \t\t ABANDON\t\t ABANDON/(NO ABANDON)%\t\t SPEED UP ");
	  System.out.println(df.format(time1)+"\t\t\t"+df.format(time2)+"\t\t\t"+(int)(100.0*(time2/time1))+"%"+"\t\t\t"+df.format(time1/time2));
	 
	  }
	  } 
	  catch (Exception ex) {
	  System.out.println("Error performing the shapelet transform"+ex);
	  ex.printStackTrace(); System.exit(0); }
	  
	  }*/
	 
}
