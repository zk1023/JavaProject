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
 * @date 2018年3月14日
 * Description
 */
public class FileUtils {
	/**
	 * 
	 * 2018年3月14日
	 * @param path
	 * @return
	 * @throws Exception
	 * Description   获取传感器个数
 	 */
	public static int getNumber_Sensor(String path) throws Exception{
		File file = new File(path) ;
		if(file.exists()&&file.isFile()){
			CSVReader reader=new CSVReader(new FileReader(path));
			List<String[]> list = reader.readAll();
			reader.close();
			int number = list.get(0).length - 1 ;
			return  number;
		}
		return 0 ;
	}
	/**
	 * 
	 * 2018年3月14日
	 * @param path
	 * @return
	 * @throws Exception
	 * Description   获取文件总行数
 	 */
	public static int getNumber_File(String path) throws Exception{
		File file = new File(path) ;
		if(file.exists()&&file.isFile()){
			CSVReader reader=new CSVReader(new FileReader(path));
			List<String[]> list = reader.readAll();
			reader.close();
			int number = list.size() ;
			return  number;
		}
		return 0 ;
	}
	/**
	 * 
	 * 2018年3月14日
	 * @param path
	 * @return
	 * @throws Exception
	 * Description   获取动作个数
 	 */
	public static int getNumber_Act(String path) throws Exception{
		File file = new File(path) ;
		if(file.exists()){
			CSVReader reader=new CSVReader(new FileReader(path));
			List<String[]> list = reader.readAll();
			reader.close();
			TreeSet<String> setKinds = new TreeSet<String>() ;
			String label ;
			for(int i = 0; i < list.size(); i ++){
				label = list.get(i)[list.get(i).length - 1];
				setKinds.add(label) ;
			}
			return  setKinds.size() ;
		}
		return 0 ;
	}
	
	/**
	 * 
	 * 2018年3月14日
	 * @param path
	 * @return
	 * @throws Exception
	 * Description   获取动作具体类别
 	 */
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
	
	/**
	 * 
	 * 2018年3月13日
	 * @param filePath
	 * Description 删除文件夹下的所有文件夹
	 */
	public static void deleteFile(String path){
		File file =  new File(path) ;
		if(file.exists()){
			String []content = file.list() ;
			if(content.length == 0){
				return ;
			}
			for(String fileName : content){
				File currentFile = new File(path + "/" + fileName) ;
				delDir(path + "/" + fileName);
				currentFile.delete() ;
			}
		}
	}
	/**
	 * 
	 * 2018年3月13日
	 * @param path 
	 * @param name
	 * Description 根据文件路径和文件名创建文件夹
	 */
	public static void createFoder(String path, String name){
		File file = new File(path + "/" + name) ;
		if(file.exists()){
			delDir(path + "/" + name);
		}else{
			file.mkdir() ;
		}
	}
	public static void createFoder(String filename){
		File file = new File(filename) ;
		if(file.exists()){
			delDir(filename);
		}else{
			file.mkdir() ;
		}
	}
	/**
	 * 
	 * 2018年3月13日
	 * @param path
	 * @param name
	 * Description 在path路径下创建n个文件夹
	 */
	public static void createFiles(String path, String name[], int n){
		for(int i = 0; i < n; i ++){
			createFoder(path, name[i]);
		}
	}
	
	/**
	 * 
	 * 2018年3月11日
	 * @param filename
	 * Description 删除文件
	 */
	public static int delFile(String filename){
        File file=new File(filename);
//        System.out.println(file.exists());
        if(file.exists()&&file.isFile()){
        	file.delete();
        	return 1 ;
        }
        return 0 ;
    }
	/**
	 * 
	 * 2018年3月12日
	 * @param filename
	 * Description 删除文件夹下所有文件
	 */
	public static void delDir(String filename){
        File file=new File(filename);
        if(file.exists()&&file.isDirectory()){
        	String[] content = file.list();//取得当前目录下所有文件
        	if(content.length == 0)
        		return ;
            for(String name : content){  
                delFile(filename + "/" + name);
        }    
       }else{
    	   createFoder(filename);
       }
	}
	/**
	 * 
	 * 2018年3月15日
	 * @param file
	 * @return
	 * Description 判断一个文件夹是否为空
	 */
	public static boolean isFileEmpty(File file) {
		if(file.isDirectory()){
			if(file.listFiles().length == 0){
				return true ;
			}
			return false ;
		}
		return false ;
	}
}
