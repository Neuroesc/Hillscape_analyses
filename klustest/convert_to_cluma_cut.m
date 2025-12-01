




for i = 1:4
    c1 = getcut(sprintf('kwiktint_%d.cut',i));

    % write an empty .cut file
    outname = sprintf('cluma.t%d.clusters',i);
    fid = fopen(outname,'w+');   
    fprintf(fid,'%d\n',c1(:));   
    fclose(fid);
end

































