/**
 * 
 */
package com.lw.process;

import java.util.ArrayList;
import java.util.List;

import java_cup.internal_error;


/**  
 * <p>Description: </p>  
 */
/**
 * @author Zhangkai
 * @date 2018年3月11日
 * Description 训练处理两条龙
 */
public class Main {
	//flag=0 代表训练（此时只会生成shapelet但是不会生成训练矩阵）	1 代表生成训练矩阵 并且 开始测试，评估模型准确度
	public static int flag =1 ;
	//从数据中抽出一部分进行训练或测试，每个传感器类别中抽取part个数据
	public static int part = 30 ;
	//时间序列长度
	public static int timeSeries = 900 ;
	//训练或者测试数据路径
	public static String path = "" ;
	//数据集路径
	public static String dataSet_path = "" ;
	//训练集测试集路径
	//数据集路径
	public static String filePath_source_train = "DataSet/Source_train" ;
	public static String filePath_source_test = "DataSet/Source_test" ;
	public static String filePath_train = "DataSet/Train" ;
	public static String filePath_test = "DataSet/Test" ;
    public static String filePath_practice_train = "DataSet_part/Train" ;
	public static String filePath_practice_test = "DataSet_part/Test" ;
	
	//训练矩阵
	public static String fileName_train = "ResultMatrix/matrixTrain.csv" ;
	//测试矩阵
	public static String fileName_test = "ResultMatrix/matrixTest.csv" ;		
	//按照动作分类后的文件位置
	public static String filePath_act = "data/wrist" ;
	//是否对训练集数据进行选择  1代表选择部分数据进行训练或测试  0代表不对数据进行选择
	public static int chooseData = 1;
	//传感器种类个数
	public static int number_sensor = 0 ;
	//动作种类个数
	public static int number_act = 0 ;
	//标签类别
	public static List<String> label = new ArrayList<String>() ;
	
	public static void main(String[]args) throws Exception {
		for(int i = 0;i<2;i++){
			if(i==0){
				flag = 0;
				part = 20;
			}else if(i==1){
				flag = 1;
				part = 50;
			}
		config();
		//把数据按动作种类整理
		getSeries.SplitClass() ;
		//按照传感器种类整理，然后转化成时间序列的数据
		getSeries.SplitSensor() ;
		//根据传感器对数据进行划分
		SensorColumn.senColumn() ;
		//给数据添加属性，使weka能够处理
		wrist_i_column_addAttribute.addAttribute();
//		选择数据和生成shapelet
		UnivariateShapelet.chooseData();
		if(flag == 1){
			getInstances.getMatrix(flag);
			AttributeSelectionTest.showResult() ;
		}
	}
	}
	
	public static void config() throws Exception{
		//默认原始训练集路径
//		path = "D:/code/Java/experiment/DataSet/subject101.csv" ;
//		path = "D:/code/Java/experiment/subject101.csv" ;
	//	path = "D:/code/Java/experiment/DataSet/Participant_1_wrist.csv" ;
//		path = "D:/code/Java/experiment/DataSet/wrist01200.csv" ;
		path = "leftsensor17_train.csv";
		//获取活动具体类别
		label = FileUtils.getLabels(path) ;
		if(label == null || label.isEmpty() || label.size() == 0 || label.size() > 1000){
			System.out.println("something is wrong 其实要么没标签，要么是文件有问题！");
			if(label != null){
				System.out.println("label个数为：     " +label.size());
			}
			return ;
		}
		//获取动作个数
		number_act  = label.size() ;
		System.out.println("动作种类：");
		for(int count = 0; count < number_act; count ++){
			System.out.print(label.get(count) + "  ");
			if(count == number_act - 1){
				System.out.println();
			}
		}
		//获取传感器个数
		number_sensor = FileUtils.getNumber_Sensor(path) ;
		System.out.println("传感器种类：\n" + number_sensor + "\n");
		if(number_sensor == 0){
			System.out.println("something is wrong 竟然没有标签，文件错误！");
			return ;
		}
		System.out.println("数据预处理开始=============================================");
		if(flag == 1){
			getInstances.getMatrix(0);
//			path = "D:/code/Java/experiment/DataSet/subject102.csv" ;
//			path = "D:/code/Java/experiment/DataSet/subject102_test_filter.csv" ;
		//	path = "D:/code/Java/experiment/subject102.csv" ;
		//	path = "D:/code/Java/experiment/DataSet/Participant_2_wrist.csv" ;
//			path = "D:/code/Java/experiment/DataSet/wrist12002400.csv" ;
			  path = "leftsensor17_test.csv" ;
		}
		System.out.println("获取数据源：" + path);
	}
}
