package com.lw.process;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class wrist_i_column_addAttribute {

	public static void addAttribute() throws IOException {
		String path = "" ;
		String oprate = "" ;
		if(Main.flag == 0){
			path = Main.filePath_source_train ;
			oprate = Main.filePath_train ;
		}else {
			path = Main.filePath_source_test ;
			oprate = Main.filePath_test ;
		}
		FileUtils.delDir(oprate);
		//修改传感器个数
		for(int i = 0; i < Main.number_sensor; i ++){
			String fileReaderName = path + "/" + "sensor_"+i+".csv" ;
			String fileWriterName = oprate + "/" + "sensor_"+i+".csv" ;
			addAttribute(fileReaderName, fileWriterName);
		}
		System.out.println("数据格式化完毕！");
	}

	private static void addAttribute(String fileReaderName,String FileWriterName) throws IOException {
		FileWriter writer = new FileWriter(FileWriterName);
        FileReader reader = new FileReader(fileReaderName);
		BufferedWriter bw = new BufferedWriter(writer);
        BufferedReader br = new BufferedReader(reader);
		bw.write("@relation" + " " + "test" + "\r\n");
		//与时间序列长度有关
		for (int i = 1; i < Main.timeSeries + 1; i++) {
			bw.write("@attribute" + " " + "attribute" + i + " " + "numeric" + "\r\n");
		}
		bw.write("@attribute" + " " + "class" + " " + "{");
		//与
		for (int i = 0; i < Main.number_act; i++) {
			if (i == (Main.number_act - 1)) {
				bw.write(Main.label.get(i) + "}" + "\r\n");
				break;
			}
			bw.write(Main.label.get(i) + ",");
		}
		bw.write("@data" + "\r\n");
		String str = "";
		while(( str = br.readLine())!=null){
			bw.write(str);
			bw.write("\n");
		}
		bw.close();
		br.close();
	}
}
