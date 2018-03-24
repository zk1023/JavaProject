package com.lw.process;
import weka.attributeSelection.*;
import weka.core.*;
import weka.classifiers.*;
import weka.classifiers.meta.*;
import weka.classifiers.trees.*;
import weka.filters.*;

import java.io.FileReader;
import java.io.FileWriter;
import java.util.*;

import com.opencsv.CSVWriter;
 
/**
 * performs attribute selection using CfsSubsetEval and GreedyStepwise
 * (backwards) and trains J48 with that. Needs 3.5.5 or higher to compile.
 *
 * @author FracPete (fracpete at waikato dot ac dot nz)
 */
public class AttributeSelectionTest {
 
  /**
   * uses the meta-classifier
   */
  protected static void useClassifier(Instances data) throws Exception {
    System.out.println("n1. Meta-classfier");
    AttributeSelectedClassifier classifier = new AttributeSelectedClassifier();
    CfsSubsetEval eval = new CfsSubsetEval();
    GreedyStepwise search = new GreedyStepwise();
    search.setSearchBackwards(true);
    RandomForest base = new RandomForest();
    classifier.setClassifier(base);
    classifier.setEvaluator(eval);
    classifier.setSearch(search);
    Evaluation evaluation = new Evaluation(data);
    evaluation.crossValidateModel(classifier, data, 10, new Random(1));
    System.out.println(evaluation.toSummaryString());
  }
 
  /**
   * uses the filter
   */
  protected static void useFilter(Instances data) throws Exception {
    System.out.println("n2. Filter");
    weka.filters.supervised.attribute.AttributeSelection filter = new weka.filters.supervised.attribute.AttributeSelection();
    CfsSubsetEval eval = new CfsSubsetEval();
    GreedyStepwise search = new GreedyStepwise();
    search.setSearchBackwards(true);
    filter.setEvaluator(eval);
    filter.setSearch(search);
    filter.setInputFormat(data);
    Instances newData = Filter.useFilter(data, filter);
   
    CSVWriter writer = new CSVWriter(new FileWriter("filterShapelet.arff"));
    for(int i = 0;i<newData.numInstances();i++){
    	writer.writeNext(newData.instance(i).toString().split(","));
    }
    writer.close();
//    RandomForest classifier = new RandomForest();
    System.out.println("***************************");
//    System.out.println(newData);
    Evaluation evaluation = new Evaluation(newData);
//    evaluation.evaluateModel(classifier, ) ;
//    evaluation.crossValidateModel(classifier, newData, 10, new Random(1));
    System.out.println("=========================");
    System.out.println(evaluation.toSummaryString());
  }
 
  /**
   * uses the low level approach
   */
  protected static void useLowLevel(Instances data) throws Exception {
    System.out.println("n3. Low-level");
    AttributeSelection attsel = new AttributeSelection();
    CfsSubsetEval eval = new CfsSubsetEval();
    GreedyStepwise search = new GreedyStepwise();
    search.setSearchBackwards(true);
    attsel.setEvaluator(eval);
    attsel.setSearch(search);
    attsel.SelectAttributes(data);
    int[] indices = attsel.selectedAttributes();
    System.out.println("selected attribute indices (starting with 0):n" + Utils.arrayToString(indices));
  }
  
  protected static void useFilter(Instances trainData,Instances testData) throws Exception {
	    System.out.println("n2. Filter");
	    weka.filters.supervised.attribute.AttributeSelection filter = new weka.filters.supervised.attribute.AttributeSelection();
	    CfsSubsetEval eval = new CfsSubsetEval();
	    GreedyStepwise search = new GreedyStepwise();
	    search.setSearchBackwards(true);
	    filter.setEvaluator(eval);
	    filter.setSearch(search);
	    filter.setInputFormat(trainData);
	    Instances newData = Filter.useFilter(trainData, filter);
	   
	   /* CSVWriter writer = new CSVWriter(new FileWriter("filterShapelet.arff"));
	    for(int i = 0;i<newData.numInstances();i++){
	    	writer.writeNext(newData.instance(i).toString().split(","));
	    }*/
	   // writer.close();
	    System.out.println("***************************");
	    RandomForest classifier = new RandomForest();
	    classifier.buildClassifier(newData);
	    Evaluation evaluation = new Evaluation(newData);
	    evaluation.evaluateModel(classifier,testData);
	  
	    // evaluation.crossValidateModel(classifier, instances_test, 10, new Random(1));
	    System.out.println(evaluation.toSummaryString());
	    System.out.println(evaluation.toMatrixString());
	   // for(int i = 0;i<testData.numInstances();i++)
	   // System.out.println(i+":"+classifier.classifyInstance(testData.instance(i))+"   "+"true value = "+testData.instance(i).classValue());
	    // 2. filter
	   
	 
	  }

  
 
  /**
   * takes a dataset as first argument
   *
   * @param args        the commandline arguments
   * @throws Exception  if something goes wrong
   */
  public static void showResult() throws Exception {
	// load data
	    System.out.println("n0. Loading data");
	   // DataSource source = new DataSource("matrix.arff");
	    //Instances data = source.getDataSet();
	    String train = Main.fileName_train; ;
		  String test = Main.fileName_test ;
	  FileReader  reader = new FileReader(test);
	//  FileReader reader_test = new FileReader(test);
	  Instances instances = new Instances(reader);
	//  Instances instances_test = new Instances(reader_test);
	    if (instances.classIndex() == -1){
	    	System.out.println("data.numAttributes() = "+instances.numAttributes());
	    	instances.setClassIndex(instances.numAttributes() - 1);
	    //	instances_test.setClassIndex(instances_test.numAttributes() - 1);
	    }
	    
	 
	    // 1. meta-classifier
	   // useClassifier(data);
	    RandomForest classifier = new RandomForest();
	  //  classifier.buildClassifier(instances);
	    Evaluation evaluation = new Evaluation(instances);
	  //  evaluation.evaluateModel(classifier,instances_test);
	  //  evaluation.crossValidateModel(classifier, instances, 10,);
	     evaluation.crossValidateModel(classifier, instances, 10, new Random(1));
	    System.out.println(evaluation.toSummaryString());
	    System.out.println(evaluation.toMatrixString());
	//    System.out.println(evaluation.toSummaryString());
//	    for(int i = 0;i<instances_test.numInstances();i++)
//	    	System.out.println(i+":"+classifier.classifyInstance(instances_test.instance(i))+"   "+"true value = "+instances_test.instance(i).classValue());
	    // 2. filter
	   // useFilter(instances,instances_test);
	 
	    // 3. low-level
	  //  useLowLevel(data);  }
  }
  
  
  public static void main(String []args) throws Exception{
	  showResult();
  }
}