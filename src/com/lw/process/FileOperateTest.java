/**
 * 
 */
package com.lw.process;

import java.io.File;
import java.io.FileReader;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;
import com.opencsv.CSVReader;


/**  
 * <p>Description: </p>  
 */
/**
 * @author Zhangkai
 * @date 2018年3月13日
 * Description
 */
public class FileOperateTest {
	public static void main(String []args) throws Exception {
//		String sourcePath = "D:/code/Java/experiment/DataSet/person1_wrist.csv" ; 
//		String filePath = "data/wrist" ;
		
		for(int i = 0; i < 3; i ++){
			ArrayList<String> arr = new ArrayList<>();
			for(int j = 0; j < 5; j ++){
				arr.add("aaa") ;
			}
			System.out.println(arr.size());
		}	
	}
	
	public static List<String> getLabels(String path) throws Exception{
		File file = new File(path) ;
		if(file.exists()){
			CSVReader reader=new CSVReader(new FileReader(path));
			List<String[]> list = reader.readAll();
			reader.close();
			TreeSet<String> setKinds = new TreeSet<String>() ;
			String label ;
			for(int i = 0; i < list.size(); i ++){
				label = list.get(i)[list.get(i).length - 1] ;
				setKinds.add(label) ;
			}
			List<String> labels = new ArrayList<String>() ;
			for(String str : setKinds){
				labels.add(str) ;
			}
			return labels ;
		}
		return null ;
	}
}
