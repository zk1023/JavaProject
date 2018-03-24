package com.lw.process;

import java.util.ArrayList;

public class MyInstance {
	String MyInstanceClass;
	ArrayList<String> content;
	public MyInstance(ArrayList<String> content ,String MyInstanceClass){
		this.MyInstanceClass = MyInstanceClass;
		this.content = content;
	}
	
	public ArrayList<String> getContent(){
		return content;
	}
	public String getInstanceClass(){
		return MyInstanceClass;
	}
	

}
