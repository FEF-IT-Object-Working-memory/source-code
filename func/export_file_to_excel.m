function export_file_to_excel(x,name,type)
% csvwrite(strcat('C:\Users\Ehsan\Desktop\PDF\',name,'.csv'),x,2,1)
file_name=strcat(cd,'\SourceFile.xlsx');
switch type
      case 0  %% bar fig1b
       a=[]; a={'All', 'M1' ,'M2'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3);           
       end

xlswrite(file_name,a,name,'A1')
    case 1  %% scatter fig S1
       a=[]; a={'monkey','Wr','Cr'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3);           
       end

xlswrite(file_name,a,name,'A1')
    case 2  %% scatter saccad loc figS1
       a=[]; a={'monkey','X','Y'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3);           
       end

xlswrite(file_name,a,name,'A1')
 case 3  %% nice plot eye figS1
       a=[]; a={'time','X in','Y in','X Out','Y Out','X wr','Y wr'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3); 
           a{i+1,4}=x(i,4);
           a{i+1,5}=x(i,5);
           a{i+1,6}=x(i,6);
           a{i+1,7}=x(i,7);
          
       end

xlswrite(file_name,a,name,'A1')
case 4  %% microsaccade  eye Fig S1
       a=[]; a={'time','micsac cr','micsac wr'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3);           
          
       end
xlswrite(file_name,a,name,'A1')
case 5  %% niceplots fig 1c firing rate
     a=[]; a={'time','FEF in','FEF OUT','IT Pref','IT Inter','IT NPref','FEF wr','IT wr'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3); 
           a{i+1,4}=x(i,4);
           a{i+1,5}=x(i,5);
           a{i+1,6}=x(i,6);
           a{i+1,7}=x(i,7);
           a{i+1,8}=x(i,8);

       end
xlswrite(file_name,a,name,'A1')
case 6  %% maps
   a=[]; a={'time','PPL Cr','PLL Wr'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3);           
          
       end
xlswrite(file_name,a,name,'A1')
case 7  %% scatter fig S1
       a=[]; a={'monkey','Out','In'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3);           
       end

xlswrite(file_name,a,name,'A1')
case 8  %% scatter fig S3
       a=[]; a={'Monkey','Visual','Delay'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3);           
       end

xlswrite(file_name,a,name,'A1')   
case 9  %% scatter fig S2b
       a=[]; a={'Center','Count%'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
       end

xlswrite(file_name,a,name,'A1') 
case 10  %% scatter fig S2d
       a=[]; a={'Neuron1','Neuron2'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
       end

xlswrite(file_name,a,name,'A1') 

case 11  %% scatter fig S2d
       a=[]; a={'monkey','AUC FEF','Diff PPL'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3);

       end

xlswrite(file_name,a,name,'A1') 
case 12  %% niceplots fig 1c firing rate
     a=[]; a={'baseline Cr','baseline Wr','Visual Cr','Visual Wr','Delay Cr','Delay Wr'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3); 
           a{i+1,4}=x(i,4);
           a{i+1,5}=x(i,5);
           a{i+1,6}=x(i,6);

       end
xlswrite(file_name,a,name,'A1')
case 13  %% maps
   a=[]; a={'Frequency',' Cr',' Wr'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3);           
          
       end
xlswrite(file_name,a,name,'D2')
case 14  %% maps
   a=[]; a={'Trial',' FEF',' IT'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3);           
          
       end
xlswrite(file_name,a,name,'A1')
case 15  %% maps
   a=[]; a={'time','High PPL ','Low PLL'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3);           
          
       end
xlswrite(file_name,a,name,'A1')
case 16  %% maps
   a=[]; a={'time','In ','Out'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3);           
          
       end
xlswrite(file_name,a,name,'A1')
case 17  %% maps
   a=[]; a={'Neurons','IT selectivity ','High-Low AUC'};
       for i=1:size(x,1)
           a{i+1,1}=x(i,1);
           a{i+1,2}=x(i,2);
           a{i+1,3}=x(i,3);           
          
       end
xlswrite(file_name,a,name,'A1')
end