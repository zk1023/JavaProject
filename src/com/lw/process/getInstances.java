package com.lw.process;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

import java_cup.internal_error;
import weka.core.Instances;
import weka.filters.timeseries.shapelet_transforms.ShapeletTransform;
public class getInstances {

	public static ShapeletTransform st = new ShapeletTransform();
	public static void getMatrix(int flag) throws Exception {
		if(flag == 0){
			System.out.println("训练矩阵生成ing");
		}else{
			System.out.println("测试矩阵生成ing");
		}
		String testPath = Main.filePath_practice_test ;
		String trainPath = Main.filePath_practice_train ;
		if(Main.chooseData == 0){
		//	testPath = Main.filePath_train ;
		//	trainPath = Main.filePath_test ;
			
			testPath = Main.filePath_test ;
			trainPath = Main.filePath_train ;
		}
		String path = "" ;
		String result = "" ;
		if(flag == 0){
			path = trainPath ;
			result = Main.fileName_train ;
		}else if(flag == 1){
			path = testPath ;
			result = Main.fileName_test ;
		}else{
			return ;
		}
		int shapelet_number = FileUtils.getNumber_File("shapelet.csv") ;
		FileUtils.delFile(result);
		FileWriter writer = new FileWriter(result);  
        BufferedWriter bw = new BufferedWriter(writer); 
        bw.write("@relation" + " " + "test" + "\r\n");     
        for (int i = 1; i < shapelet_number + 1; i++) {  
            bw.write("@attribute" + " " + "attribute" + i + " " + "numeric" + "\r\n");  
        }  
        bw.write("@attribute"+" "+"class"+" "+"{"); 
        //添加属性 与动作种类有关
        for(int i =0;i < Main.number_act;i++){  
            if(i==(Main.number_act - 1)){  
                bw.write(Main.label.get(i)+"}"+"\r\n");  
                break;  
            }  
            bw.write(Main.label.get(i)+",");  
        }  
        bw.write("@data" + "\r\n"); 
		//传感器个数
		File files[] = new File[Main.number_sensor] ;
		FileReader fileReader[] = new FileReader[Main.number_sensor] ;
		Instances instances[] = new Instances[Main.number_sensor] ;
		for(int i = 0; i < Main.number_sensor; i ++){
			files[i] = new File(path + "/" + "sensor_" + i + ".csv") ;
//			files[i] = new File(path + "/" + "wrist_" + i + "_column.csv") ;
			fileReader[i] = new FileReader(files[i]) ;
			instances[i] = new Instances(fileReader[i]) ;
		}
		String lineStr = "";
		String str = null;
        String line = new String();
		for(int i = 0; i <instances[0].numInstances(); i ++){
			ArrayList<String> arr = new ArrayList<>();
			for(int j = 0; j < Main.number_sensor; j ++){
				lineStr = instances[j].instance(i).toString() ;
				arr.add(lineStr) ;
			}
			if(lineStr.isEmpty()){
				System.out.println("获取实例失败");
				bw.close();
				return ;
			}
		//	System.out.println("============================");
		//	System.out.println(lineStr);
			MyInstance myInstance = new MyInstance(arr, lineStr.split(",")[lineStr.split(",").length - 1]);
			BufferedReader bufferedReader = new BufferedReader(new FileReader(new File("shapelet.csv")));
			while ((str = bufferedReader.readLine()) != null) {
				str=str.replace("\"", "");
				for(int k = 0; k < Main.number_sensor; k ++){
					if (str.split(",")[str.split(",").length - 2].equals(Integer.toString(k))) {
					//	line = line +(new DTW().getDistance(myInstance.content.get(k), str)+"")+",";
						line = line+st.subsequenceDistance2(myInstance.content.get(k), str)+",";
					//System.out.println("DTW = "+new DTW().getDistance(myInstance.content.get(k), str));
					}
				}
			}
			line =line + myInstance.MyInstanceClass+"";;
			bw.write(line+"\n");
			line = "";
			bufferedReader.close();			
		}
		bw.close();
		if(flag == 0){
			System.out.println("训练矩阵生成完毕！");
		}else{
			System.out.println("测试矩阵生成完毕！");
		}
	}
}
