import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;

public class aaaa {

	public static void main(String[] args) throws IOException{
		CSVReader reader = new CSVReader(new FileReader("smartphoneatwrist.csv"));
	    List<String[]> list = reader.readAll();
	    List<String[]> subList = getData(list,2400,4800);
	    reader.close();
	    CSVWriter writer = new CSVWriter(new FileWriter("wrist.csv"));
	    writeFile(writer,subList);
	}
	private static void writeFile(CSVWriter writer, List<String[]> subList) throws IOException {
		// TODO Auto-generated method stub
		for(int i = 0;i<subList.size();i++){
			writer.writeNext(subList.get(i));
		}
		writer.flush();
		writer.close();
	}
	
	public static List<String[]> getData(List<String[]> list,int start,int length){
		List<String[]> subList = new ArrayList<>();
		int count = 0;
		for(int i = start;i<length;i++){
			subList.add(list.get(i));
		}
		for(int i = 9000+start;i<9000+length;i++){
			subList.add(list.get(i));
		}
		for(int i = 18000+start;i<18000+length;i++){
			subList.add(list.get(i));
		}
		for(int i = 27000+start;i<27000+length;i++){
			subList.add(list.get(i));
		}
		for(int i = 36000+start;i<36000+length;i++){
			subList.add(list.get(i));
		}
		for(int i = 45000+start;i<45000+length;i++){
			subList.add(list.get(i));
		}
		for(int i = 54000+start;i<54000+length;i++){
			subList.add(list.get(i));
		}
		for(int i = 63000+start;i<63000+length;i++){
			subList.add(list.get(i));
		}
		for(int i = 72000+start;i<72000+length;i++){
			subList.add(list.get(i));
		}
		for(int i = 81000+start;i<81000+length;i++){
			subList.add(list.get(i));
		}
		for(int i = 90000+start;i<90000+length;i++){
			subList.add(list.get(i));
		}
		for(int i = 99000+start;i<99000+length;i++){
			subList.add(list.get(i));
		}
		
		return subList;
		
	}
}
