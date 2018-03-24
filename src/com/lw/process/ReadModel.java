package com.lw.process;

import java.io.FileReader;

import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.core.Instances;

public class ReadModel {
	public static void main(String[] args) throws Exception{
		Instances train1=null,test1=null;
		FileReader r1;
		 Evaluation evaluation_RandomForest = null ;
		r1= new FileReader("classtest.csv");       			
			test1 = new Instances(r1); 
			test1.setClassIndex(test1.numAttributes()-1); 
		Classifier classifier8 = (Classifier) weka.core.SerializationHelper.read("LibSVM.model");  
		//evaluation_RandomForest.evaluateModel(classifier8, test1);
		for(int i = 0;i<test1.numInstances();i++)
		{
			System.out.println(classifier8.classifyInstance(test1.instance(i)));
		}
	}

}
