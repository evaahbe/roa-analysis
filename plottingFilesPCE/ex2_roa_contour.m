% this script only works for 2D stochastic systems!

function ex2_roa_contour(filenames)

    close all

    cd(filenames.resultsdirectory)
    k = 0;
    for i = 1:length(filenames.filenamestoplot)
        file_list = dir('*.mat'); %load the file names into a list
        for j = 1:length(file_list)
            cur_file_name = file_list(j).name;
            if (strfind(cur_file_name,filenames.filenamestoplot) == 1)
                k = k+1;
                filestoplot_list{k} = cur_file_name;
            end
        end
    end
    cd ..
    

    for jj = 1:length(filestoplot_list)
        
        load(filestoplot_list{jj})

        colors = getCustomColors();

        domain = [ -2*sqrt(1/Q0(1,1)) 2*sqrt(1/Q0(1,1)) -2*sqrt(1/Q0(2,2)) 2*sqrt(1/Q0(2,2))]; %this is a first rough estimate, adjust if needed.
        Nx = 200;
        Ny = 200;
        xgp = linspace(domain(1),domain(2),Nx);
        ygp = linspace(domain(3),domain(4),Ny);
        [xg,yg] = meshgrid(xgp,ygp);
        
        lypgrid = zeros(1,length(xg(:)));

        for j=1:length(xg(:))
            x01 = xg(j);
            x02 = yg(j);
            lypgrid(j) = [x01;x02;x01^2;x01*x02;x02^2]'* Q0* [x01;x02;x01^2;x01*x02;x02^2];
        end

        lypgrid = reshape(lypgrid,size(xg));
        vlyp = [1 1]; %if alpha=1;
        contourmatrix = contourc(xgp,ygp,lypgrid,vlyp);
        contourmatrix = contourmatrix(:,2:end);
        domlim = max(abs(domain));
        contourmatrix = (abs(contourmatrix)<=domlim).*contourmatrix + ((contourmatrix)>domlim).*domlim  -((contourmatrix)<-domlim).*domlim;
            
        plot(contourmatrix(1,:)+xiep0(1),contourmatrix(2,:)+xiep0(1),'LineWidth',2,'Color',colors.mycolors(jj,:))
 
        save(filestoplot_list{jj}, '-append','contourmatrix');
    end   
    


end




    
    

    
    
