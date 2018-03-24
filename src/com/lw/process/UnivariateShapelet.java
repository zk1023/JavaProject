package com.lw.process;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import com.lw.shapelet.QualityMeasures;
import com.lw.shapelet.Shapelet;
import com.opencsv.CSVWriter;

import weka.core.Instances;
import weka.filters.timeseries.shapelet_transforms.ClusteredShapeletTransform;
import weka.filters.timeseries.shapelet_transforms.ShapeletTransform;

public class UnivariateShapelet {

	public static ShapeletTransform st;
	public static ArrayList<Shapelet> clusteredShapelets;
	public static void chooseData() throws Exception{
		String path_from = "" ;
		String path_to = "" ;
		if(Main.flag == 0){
			path_from = Main.filePath_train ;
			path_to = Main.filePath_practice_train ;
		}else if(Main.flag == 1){
			path_from = Main.filePath_test ;
			path_to = Main.filePath_practice_test ;
		}
		if(Main.chooseData == 1){
			System.out.println(path_from + "---->" + path_to);
			getData(path_from, path_to);
		}
		if(Main.flag == 0){
			generateShapelet();
		}
	}
	private static void writeShapelet(ArrayList<Shapelet> clusteredShapelets, int a) throws IOException {
		// TODO Auto-generated method stub
		// FileWriter writer = new FileWriter("shapelet.csv",true);
		CSVWriter writer = new CSVWriter(new FileWriter("shapelet.csv", true));
		for (int i = 0; i < clusteredShapelets.size(); i++) {
			double[] shapelet = clusteredShapelets.get(i).getContent();
			String[] shapelet_string = new String[shapelet.length + 2];
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
	public static Instances basicTransformExample(Instances train) {
		st = new ShapeletTransform();
		System.out.println("classsssssssssssss = "+train.instance(0).classValue());
		// Let m=train.numAttributes()-1 (series length)
		// Let n= train.numInstances() (number of series)
		int nosShapelets = (train.numAttributes() - 1) * train.numInstances() / 5;
		if (nosShapelets < ShapeletTransform.DEFAULT_NUMSHAPELETS)
			nosShapelets = ShapeletTransform.DEFAULT_NUMSHAPELETS;
		st.setNumberOfShapelets(nosShapelets);
		int minLength =(train.numAttributes() - 1) / 8;
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
		//System.out.println("----------"+st.shapelets.size());
		System.out.println(" Clustering down to " + nosShapelets + " Shapelets");
		System.out.println(" From " + st.getNumberOfShapelets() + " Shapelets");
		// ShapeletTransform st = new ShapeletTransform();
		try {
		//	System.out.println(train);
		//	System.out.println("------------"+cst.allShapelets.size());
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
	/**
	 * 由于数据量庞大，提取出部分数据
	 * @throws Exception 
	 */
	public static void getData2(String path_from, String path_to) throws Exception{
		FileUtils.delDir(path_to);
		int num = FileUtils.getNumber_File(Main.path) / Main.number_act / Main.timeSeries ;
		//与传感器个数有关
		for(int  i = 0; i< Main.number_sensor; i ++){
			BufferedReader reader = new BufferedReader(new FileReader(new File(path_from + "/sensor_"+i+".csv")));
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(path_to + "/sensor_"+i+".csv")));
			//与时间序列长度有关
			for(int j = 0; j < Main.timeSeries + 3; j ++){
			     String atributes = reader.readLine();
			     writer.write(atributes);
			     writer.write("\n");
			}
			String str = "";
			//与动作种类有关
			for(int a = 0 ;a < Main.number_act; a++){
				for(int count = 0; count < Main.part; count ++){
					str = reader.readLine();
					writer.write(str);
					writer.write("\n");
				}
				for(int count = 0; count < num - Main.part; count++){
					reader.readLine();
				}
			} 
			writer.flush();
			reader.close();
			writer.close();
		}
		System.out.println("数据选择完毕！");
	}
	
	
	public static void getData(String path_from, String path_to) throws Exception {
		FileUtils.delDir(path_to);
		
		//与传感器个数有关
		for(int  i = 0; i< Main.number_sensor; i ++){
			BufferedReader reader = new BufferedReader(new FileReader(new File(path_from + "/sensor_"+i+".csv")));
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(path_to + "/sensor_"+i+".csv")));
			String str = "";
			int mark = 0 ;
			int count[] = new int[Main.number_act] ;
			for(int j = 0; j < Main.number_act; j ++){
				count[j] = 0;
			}
			while ((str = reader.readLine()) != null){
				str = str.replace("\"", "") ;
				String[] strs = str.split(",");
				if(mark < Main.timeSeries + 3){
				     writer.write(str);
				     writer.write("\n");
				     mark ++ ;
				}else {
//					System.out.println("part"+Main.part+"label  " + strs[Main.timeSeries]);
					for(int k = 0; k < Main.number_act; k ++){
						if (strs[Main.timeSeries].equals(Main.label.get(k)) && (count[k] < Main.part)) {
							count[k] ++;
							writer.write(str + "\n") ;
						}
					}
				}
			}
			writer.close();
			reader.close();
		}
		System.out.println("数据选择完毕！");
	}
	
	
	
	public static void generateShapelet() throws Exception{
		String path = Main.filePath_practice_train ;
		if(Main.chooseData == 0){
			path = Main.filePath_train ;
		}
		FileUtils.delFile("shapelet.csv") ;
		//传感器个数
		for(int i = 0; i < Main.number_sensor; i ++){
			System.out.println("生成shapeleting:============================= " + i);
			//根据训练集中的文件生成shapelet
			FileReader reader = new FileReader(new File(path+ "/" + "sensor_"+i+".csv"));
//			FileReader reader = new FileReader(new File(path+ "/" + "wrist_"+i+"_column.csv"));
			Instances train = new Instances(reader);
			train.setClassIndex(train.numAttributes()-1);
			basicTransformExample(train);
			clusteredShapeletTransformExample(train);
			writeShapelet(clusteredShapelets, i);
		}
		System.out.println("Shapelet生成完毕！");
	}
	
	/**
	 * Description: 清除shapelet格式
	 * 
	 */
	public static void clearFormat(){
		
	}
}
