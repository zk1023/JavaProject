package com.lw.process;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;

public class SensorColumn {
	public static void senColumn() throws IOException {
		File file = new File(Main.filePath_act);
		if(file.listFiles().length == 0){
			System.out.println(Main.filePath_act + " 空文件夹！ ");
			return ;
		}
		String path = "" ;
		if(Main.flag == 0){
			path = Main.filePath_source_train ;
		}else {
			path = Main.filePath_source_test ;
		}
		FileUtils.delDir(path);
		String file_name ;
		File[] folder = file.listFiles();
		if(FileUtils.isFileEmpty(folder[0])){
			System.out.println(folder[0].getName() + " 空文件夹！ ");
			return ;
		}
		//传感器个数
		for (int i = 0; i < Main.number_act; i++) {
			for (int j = 0; j < Main.number_sensor; j ++) {
				file_name = "sensor_" + j +".csv" ;
				FileWriter fileWr = new FileWriter(path + "/" + file_name, true) ;
				CSVWriter csvWr = new CSVWriter(fileWr) ;
				write(csvWr, Main.filePath_act + "/" + Main.label.get(i) + "/" + file_name, Main.label.get(i));
			}
		}
		System.out.println("根据传感器类别对数据划分完毕！");
	}

	private static void write(CSVWriter writer, String filename, String class_value) throws IOException {
		CSVReader reader = new CSVReader(new FileReader(new File(filename)));
		List<String[]> list = reader.readAll();
		String[] list2 = new String[list.get(0).length + 1];
		for (int i = 0; i < list.size(); i++) {
			for (int j = 0; j < list.get(i).length; j++) {
				list2[j] = list.get(i)[j];
			}
			list2[list2.length - 1] = class_value;
			writer.writeNext(list2);
		}
		reader.close();
		writer.close();
	}
}
