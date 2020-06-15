% this script only works for 2D stochastic systems!

function ex1_plot(filenames)

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
    
    width = 40;
    height = 28;
    roaplot = figure('Units','centimeters',...
    'Position',[1 0.7 width height],...
    'PaperPositionMode','auto');
    hold on
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
        if numsetsRE.degs.V_dU == 2
            for j=1:length(xg(:))
                x01 = xg(j);
                x02 = yg(j);
                lypgrid(j) = [x01;x02]'* Q0* [x01;x02];
            end
        elseif numsetsRE.degs.V_dU == 4
            for j=1:length(xg(:))
                x01 = xg(j);
                x02 = yg(j);
                lypgrid(j) = [x01;x02;x01^2;x01*x02;x02^2]'* Q0* [x01;x02;x01^2;x01*x02;x02^2];
            end
        end
        lypgrid = reshape(lypgrid,size(xg));
        vlyp = [1 1]; %if alpha=1;
        contourmatrix = contourc(xgp,ygp,lypgrid,vlyp);
        contourmatrix = contourmatrix(:,2:end);
        domlim = max(abs(domain));
        contourmatrix = (abs(contourmatrix)<=domlim).*contourmatrix + ((contourmatrix)>domlim).*domlim  -((contourmatrix)<-domlim).*domlim;
            
        plot(contourmatrix(1,:)+xiep0(1),contourmatrix(2,:)+xiep0(1),'LineWidth',2,'Color',colors.mycolors(jj,:))

        legendentries{jj} =  strcat('$\partial(V) =$',num2str(numsetsRE.degs.V_dU),', $\hat{\sigma}_{ii} =$ [', num2str(sys.varfix(1,1)), ',', num2str(sys.varfix(2,2)), ']');

    end   
    
    
    load VDP_LCtraj_mu07.mat
    x_lc = xtraj{2};
    plot(x_lc(:,1), x_lc(:,2),':k','LineWidth',1.5) 
    legendentries{jj+1} = 'LC $(c=0.7)$';

    load VDP_LCtraj_mu13.mat
    x_lc = xtraj{2};
    plot(x_lc(:,1), x_lc(:,2),'--k','LineWidth',1.5)
    legendentries{jj+2} = 'LC $(c=1.3)$';
    

    plot(xiep0(1), xiep0(2), 'ko','LineWidth',1.5)
    legendentries{jj+3} = '$x_{EP}$';
    
    set(gca,...
    'LineWidth', 1.5,...
    'Units','normalized',...                  
    'FontSize',20)
    title('ROA estimates','FontSize',20)
    legend(legendentries,'LineWidth',1.2,'FontSize',23,'NumColumns',1,'interpreter','latex')
    xlabel('$\bar{x}_{1_0}$','FontSize',35,'interpreter','latex');
    ylabel('$\bar{x}_{2_0}$','FontSize',35,'interpreter','latex');  

end




    
    

    
    
