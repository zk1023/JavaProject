package com.lw.process;

import weka.filters.timeseries.shapelet_transforms.ShapeletTransform;

public class DTW {

	public double getMin(double a, double b, double c) {
		double min = a;
		if (b > a)
			min = a;
		else if (c > b) {
			min = b;
		} else {
			min = c;
		}
		return min;
	}
/*
 * seqa:instances
 * seqb:shapelet
 */
	public double getDistance(String seqa, String seqb) {
		double distance = 0;
		double[] seqa_discre = new double[seqa.split(",").length-1];
		for(int i =0;i<seqa.split(",").length-1;i++){
			seqa_discre[i] = Double.parseDouble(seqa.split(",")[i]);
		}
	//	seqa_discre = zNorm(seqa_discre, false);
		int lena = seqa_discre.length;
		int lenb = seqb.split(",").length-2;
		double[][] c = new double[lena][lenb];
		for (int i = 0; i < lena; i++) {
			for (int j = 0; j < lenb; j++) {
				c[i][j] = 1;
			}
		}
		for (int i = 0; i < lena; i++) {
			for (int j = 0; j < lenb; j++) {
				double tmp = ((seqa_discre[i]) - Double.parseDouble(seqb.split(",")[j])) * ((seqa_discre[i]) - Double.parseDouble(seqb.split(",")[j]));
				if (j == 0 && i == 0)
					c[i][j] = tmp;
				else if (j > 0)
					c[i][j] = c[i][j - 1] + tmp;
				if (i > 0) {
					if (j == 0)
						c[i][j] = tmp + c[i - 1][j];
					else
						c[i][j] = tmp + getMin(c[i][j - 1], c[i - 1][j - 1], c[i - 1][j]);
				}
			}
		}
		distance = c[lena - 1][lenb - 1];
		return distance;
	}
	protected double[] zNorm(double[] input, boolean classValOn){        
        return ShapeletTransform.zNormalise(input, classValOn);
    }

}