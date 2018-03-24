package com.lw.process;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;



public class getSeries {
	/**
	 * 
	 * @throws FileNotFoundException
	 * @throws IOException
	 * Description:把数据按照动作种类进行分开
	 */
	public static void SplitClass() throws Exception {
		CSVReader reader=new CSVReader(new FileReader(Main.path));
		int number = Main.number_act ;
		FileWriter fileWriter[] = new FileWriter[number] ;
		CSVWriter csvWriter[] = new CSVWriter[number] ;
		FileUtils.deleteFile(Main.filePath_act);
		for(int i = 0; i < number; i ++){
			fileWriter[i] = new FileWriter(Main.filePath_act + "/" + Main.label.get(i) + ".csv") ;
			csvWriter[i] = new CSVWriter(fileWriter[i]) ;
		}	
		List<String[]> list = reader.readAll();
		 for(int i = 0;i<list.size();i++){
			 for(int j = 0; j < number; j ++){
				 if(list.get(i)[list.get(i).length-1].equals(Main.label.get(j))){
					 csvWriter[j].writeNext(list.get(i));
					 continue ;
				 }
			 }	 
		 }
		 reader.close();
		 for(int i = 0; i < number; i ++){
			 csvWriter[i].close();
		 }
		 System.out.println("根据动作类别对数据划分完毕！");	
	}

	/**
	 * 
	 * @throws FileNotFoundException
	 * @throws IOException
	 * Description:将同一个类中的不同传感器的不同轴分开分别存到不同的文件，目的是用来提取序列
	 */
	public static void SplitSensor() throws Exception {
		int number = Main.number_act ;
		FileReader fileReader[] = new FileReader[number] ;
		CSVReader csvReader[] = new CSVReader[number] ;
		//动作种类
		for(int i = 0; i < number; i ++){
			fileReader[i] = new FileReader(Main.filePath_act + "/" + Main.label.get(i) + ".csv") ;
			csvReader[i] = new CSVReader(fileReader[i]) ;
			FileUtils.createFoder(Main.filePath_act, Main.label.get(i));
			writeToFile(csvReader[i].readAll(), Main.label.get(i)) ;
			csvReader[i].close();
		}
		 //删除文件，与动作种类有关
		 for(int i = 0; i < number; i ++){
			 FileUtils.delFile(Main.filePath_act + "/" +  Main.label.get(i) + ".csv");			 
		 }	 
		 System.out.println("把数据转化为时间序列 ，并把每中动作分别按照不同传感器进一步分类！");
	}

	
/**
 * 
 * @param list 数据
 * @param str 写入文件名
 * @throws IOException
 * Description： 时间序列转化
 */
	public static void writeToFile(List<String[]> list, String str) throws Exception {
		int block_size = Main.timeSeries ;
		for(int i = 0;i<Main.number_sensor;i++){
			 CSVWriter writer1 = new CSVWriter(new FileWriter(Main.filePath_act + "/" + str+"/" + "sensor_"+i+".csv"));
			 //后面还需要计算时间序列条数
			 for(int j = 0;j<list.size()/block_size;j++){
				 String[] value = new String[block_size];
				 for(int k = j*block_size;k<block_size*(j+1);k++){					 
					 value[k-j*block_size]= list.get(k)[i];
				 }
				 writer1.writeNext(value);
				 writer1.flush();
			 }
			 writer1.close();
		 }
	}
}
