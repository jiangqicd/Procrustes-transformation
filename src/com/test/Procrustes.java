package ClusterFinding;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class Procrustes {
    public static Matrix deviation(double[] mat) {
        Matrix m = new Matrix(mat.length,1);
        for (int i=0;i<mat.length;i++){
                m.set(i, 0, mat[i]);
        }
        return m;
    }
    public static Matrix deviation_nochange(double[][] mat) {
        Matrix m = new Matrix(mat.length,mat[0].length);
        for (int i=0;i<mat.length;i++){
            for(int j=0;j<mat[0].length;j++){
                m.set(i, j, mat[i][j]);
            }

        }
        return m;
    }
    public static Matrix deviation(double[][] mat) {
        Matrix m = new Matrix(mat.length, mat[0].length);
        for (int j=0;j<mat[0].length;j++){
            for (int i=0;i<mat.length;i++){
                m.set(i, j, mat[i][j]);
            }
        }
        return m;
    }
    public static double[] centroids(double[][] mat){
        double []centroids=new double[mat[0].length];
        for(int i=0;i<mat[0].length;i++){
            for(int j=0;j<mat.length;j++){
                centroids[i]+=mat[j][i];
            }
        }
        System.out.print("centroid:"+"\n");
        System.out.print("[");
        for(int i=0;i<centroids.length;i++){
            centroids[i]=centroids[i]/mat.length;
            System.out.print(centroids[i]+"  ");
        }
        System.out.print("]"+"\n");
        return centroids;
    }

    public static double[] move(double[] mat1,double[]mat2){
        double[] move=new double[mat1.length];
        System.out.print("move factor t:"+"\n");
        System.out.print("[");
        for(int i=0;i<mat1.length;i++){
            move[i]=mat1[i]-mat2[i];
            System.out.print(move[i]+"  ");
        }
        System.out.print("]"+"\n");
        return move;
    }

    public static double fun_c(double[][] mat1,double[][] mat2,double[] c_id1,double[] c_id2,double[] v){
        double c=1;
        Matrix Mat_1 = deviation(mat1);
        Matrix Mat_2= deviation(mat2);
        Matrix C_id1=deviation(c_id1);
        Matrix C_id2=deviation(c_id2);
        Matrix C_id1T=C_id1.transpose();
        Matrix C_id2T=C_id2.transpose();
        Matrix V=deviation(v);
        Matrix vp1=V.times(C_id1T);
        Matrix vp2=V.times(C_id2T);
        Matrix Pre=Mat_1.minus(vp1);
        Matrix Cur=Mat_2.minus(vp2);
        double a=Pre.norm2();
        double b=Cur.norm2();
        c=a/b;
        System.out.print("scaling factor c:"+"\n");
        System.out.print(c+"\n");
        return c;
    }

    public static double[][] fun_R(double[][] mat1,double[][] mat2,double[] t,double c,double[] v){
        double [][] R=new double[mat1[0].length][mat1[0].length];
        Matrix Mat_1 = deviation(mat1);
        Matrix Mat_2= deviation(mat2);
        Matrix Mat_1T=Mat_1.transpose();
        Matrix T=deviation(t);
        Matrix TT=T.transpose();
        Matrix V=deviation(v);
        Matrix vt=V.times(TT);
        Matrix m=Mat_2.plus(vt);
        Matrix cMat1_T=Mat_1T.times(c);
        Matrix Resulat=cMat1_T.times(m);
        SingularValueDecomposition svd = Resulat.svd();
        Matrix VV=svd.getV();
        Matrix UU=svd.getU();
        Matrix UUT=UU.transpose();
        Matrix R1=VV.times(UUT);
        System.out.print("R:"+"\n");
        for(int i=0;i<mat1[0].length;i++){
            for(int j=0;j<mat1[0].length;j++){
                R[i][j]=R1.get(i,j);
                System.out.print(R[i][j]+"  ");            }
            System.out.print("\n");
        }
        return R;
    }

    public static double[][] fun_P(double[][] mat1,double[] t,double c,double[] v,double[][] R){
        double [][] P=new double[mat1.length][mat1[0].length];
        Matrix Mat_1 = deviation(mat1);
        Matrix R1=deviation_nochange(R);
        Matrix T=deviation(t);
        Matrix TT=T.transpose();
        Matrix V=deviation(v);
        Matrix vt=V.times(TT);
        Matrix m=Mat_1.plus(vt);
        m.print(mat1.length,2);
        Matrix cm=m.times(c);
        cm.print(mat1.length,2);
        R1.print(2,2);
        Matrix Resulat=cm.times(R1);
        Resulat.print(mat1.length,2);

        System.out.print("P_new:"+"\n");
        for(int i=0;i<mat1.length;i++){
            for(int j=0;j<mat1[0].length;j++){
                P[i][j]=Resulat.get(i,j);
                System.out.print(P[i][j]+",");            }
            System.out.print("\n");
        }
        return P;
    }
    public static void main (String[] args){
         double[][] matrix1 = new double[20][2];// matrix1 is  previous matrix P
        System.out.print("P_previous:"+"\n");
         for(int i=0;i<20;i++){
             for(int j=0;j<2;j++){
                 matrix1[i][j]=i+j;
                 System.out.print(matrix1[i][j]+",");
             }
             System.out.print("\n");
         }
        double[][] matrix2 = new double[20][2];//matrix1 is  current matrix P'
        System.out.print("P_current:"+"\n");
        for(int i=0;i<20;i++){
            for(int j=0;j<2;j++){
                matrix2[i][j]=(i+j)*2;
                System.out.print(matrix2[i][j]+",");
            }
            System.out.print("\n");
        }

        Matrix mat = deviation(matrix1);
        SingularValueDecomposition svd = mat.svd();//svd
        double[] lambda = svd.getSingularValues();
        Matrix U=svd.getU();
        Matrix V=svd.getV();
        Matrix S=svd.getS();

        //求质心
        double[] centroid1=centroids(matrix1);
        double[] centroid2=centroids(matrix2);

        //求偏移t
        double[] t=move(centroid1,centroid2);

        //求单位向量v

        double[] v=new double[matrix1.length];
        for(int i=0;i<matrix1.length;i++){
            v[i]=1;
        }

        //求缩放因子c

        double c=fun_c(matrix1,matrix2,centroid1,centroid2,v);

        //求旋转矩阵R

        double[][] R=fun_R(matrix1,matrix2,t,c,v);

        //求新的P矩阵
        double[][] p=fun_P(matrix2,t,c,v,R);
    }
}
