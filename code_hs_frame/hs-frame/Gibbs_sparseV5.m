function [I_syn, I_squared_sum, response_IB]= Gibbs_sparseV5(template, filters, I_syn, corrBB, c_val_list, I_squared_sum, response_IB, nonzero_corrBB_List, lower_bound_rand, upper_bound_rand)

   for j=1:template.numSelected  % go over all the selected bases     
           %disp(['sampling for ' num2str(j) '-th basis of ' num2str(template.numSelected)]);
           % draw sample for c which is relevant to B_j     
           
           %rand_num =  (rand()*32767 + 1)/ (32767+2);
           rand_num =  lower_bound_rand + (upper_bound_rand - lower_bound_rand).*rand(1,1);
           c=gibbs_v3(template.selectedLambdas, response_IB, single(c_val_list), corrBB, j, single(I_squared_sum), response_IB(j), single(rand_num), nonzero_corrBB_List{j});    
    
           %c_Matrix(iter,j)=c; % save the values of c
       
           % update I by I=I+cB; 
           I_syn=I_plus_cB(I_syn, c, filters{template.selectedFilter(j)}, template.selectedRow(j), template.selectedCol(j)); 
    
           % update I_squared  I^2
           I_squared_sum = I_squared_sum + c^2 + 2*c* response_IB(j);      
           
           % update <I, B_i> for all i
%            updateElement=nonzero_corrBB_List{j};
%            for i=1:length(updateElement) 
%                ii=updateElement(i);
%                response_IB(ii)= response_IB(ii) + c* corrBB(j, ii);        
%            end     
           
           for i=1:template.numSelected
               response_IB(i)= response_IB(i) + c* corrBB(j, i);        
           end        
           
   end