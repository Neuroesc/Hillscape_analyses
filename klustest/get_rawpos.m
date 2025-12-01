function [tdata,post] = get_rawpos(posfile)

    [post,~,~,~,targets,~,~] = Nlx2MatVT(posfile,[1 1 1 1 1 1], 1, 1, [] ); 
    post = (post(:) - post(1)) ./ 1e6;

    nsamps = numel(post);    
    tdata = zeros(nsamps,2,3);
    for i =1:nsamps
        for j=1:4 % only need to use the first 4 values           
            target_value = targets(j,i);            
            binary_target_string = dec2bin(target_value,32);

            x = bin2dec(binary_target_string(21:32));
            y = bin2dec(binary_target_string(5:16));            

            if bin2dec(binary_target_string(2)) % red
                tdata(i,:,1) = [x y];
            elseif bin2dec(binary_target_string(3)) % green
                tdata(i,:,2) = [x y];
            elseif bin2dec(binary_target_string(4)) % blue
                tdata(i,:,3) = [x y];
            end
        end
    end

    
    
    
    
    
    
%     text_headers = {'trial_date','trial_time','experimenter','comments','sw_version','pos_format'};
%     digit_headers = {'duration','num_colours','min_x','max_x','min_y','max_y','window_min_x','window_max_x','window_min_y','window_max_y','bytes_per_timestamp'};
%     units_headers = {'timebase','sample_rate'};
% 
%     fid = fopen(posfile,'r','ieee-be'); % open the set file for reading
%     line = fgets(fid); % look at the next line of the file (end line = -1)
%     headers = struct;
%     ind = 1;
%     while line ~= -1 
%         if contains('data_start',line(1:10)) % test if this is the start of the data
%             break
%         else
%             [t,r] = strtok(line);
%             % text headers
%             idx = contains(text_headers,t);
%             if any( idx )
%                 headers.( text_headers{idx} ) = r(2:end);
%             end
% 
%             % numerical headers
%             idx = contains(digit_headers,t);
%             if any( idx )
%                 headers.( digit_headers{idx} ) = str2double( r(2:end) );
%             end
% 
%             % headers with units at the end
%             idx = contains(units_headers,t);
%             if any( idx )
%                 str1 = regexprep(r,'[,;=]', ' ');
%                 str2 = regexprep(regexprep(str1,'[^- 0-9.eE(,)/]',''), ' \D* ',' ');
%                 str3 = regexprep(str2, {'\.\s','\E\s','\e\s','\s\E','\s\e'},' ');
%                 headers.( units_headers{idx} ) = str2num(str3);
%             end
% 
%             % get next line
%             line = fgets(fid);  
%             ind = ind+1;        
%         end      
%     end
% 
%     fseek(fid,ind,0);
%     ind = 1;
%     temp = [];
%     while line ~= -1 
%         d = fread(fid,1,'uint32');
%         if feof(fid)
%             break
%         end    
%         post(ind) = d;
%         for jj = 1:8
%             d = fread(fid,1,'uint16');  
%             if feof(fid)
%                 break
%             end
%             temp(ind,jj) = d;  
%         end
%         ind = ind+1;
%     end    
    
    
    
    
    
    
    
    
    
    
    
    