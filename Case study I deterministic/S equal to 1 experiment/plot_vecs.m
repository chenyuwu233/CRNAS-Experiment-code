%% Description:
%  plot_vecs is a function used to plot the boxplot and add significant bar
%  according to corresponding test
%
%  Input:
%  - vec_mat: m x n matrix that records all the m x 1 vectors that needs to
%  compare.
%  - vec_name: 1 x n cell of string about the name of vectors.
%  - Y_LIM: range of y.
%  - Y_label: label for Y axis
%  - title: Title for the plot
%  - colors: n x 3 matrix that records rgb values for n box.
%  - test_cmd: command about what test to used.
%    currently supported:
%    -- rank sum (unpaired): Wilcoxon rank sum test
%    -- signed rank (paired): Wilcoxon signed rank test
%  - scale_cmd: 0 for linear 1 for log.




function plot_vecs(vec_mat,vec_name,Y_LIM,Y_label,title,colors,test_cmd,scale_cmd)
    num_vec = size(vec_mat,2);
    switch test_cmd
        case 'rank sum'
            %% Boxplotting
            boxplot(vec_mat,'Symbol','o','OutlierSize',3,'Colors',[0,0,0])
            ylim(Y_LIM)
            ax=gca;hold on;
            ax.LineWidth=1.1;
            ax.FontSize=22;
            ax.FontName='Arial';
            ax.FontWeight = 'bold';
            ax.XTickLabel=vec_name;
            ax.Title.String=title;
            ax.Title.FontSize=27;
            ax.YLabel.String=Y_label;
            if scale_cmd == 1
                ax.YScale = 'log';
            end

            %% Width of Line
            lineObj=findobj(gca,'Type','Line');
            for i=1:length(lineObj)
                lineObj(i).LineWidth=1;
                lineObj(i).MarkerFaceColor=[1,1,1].*.3;
                lineObj(i).MarkerEdgeColor=[1,1,1].*.3;
            end

            %% Coloring the box
            if ~isempty(colors)
                boxObj=findobj(gca,'Tag','Box');
                for i=1:length(boxObj)
                    patch(boxObj(i).XData,boxObj(i).YData,colors(i,:),'FaceAlpha',0.5,'LineWidth',1.1);
                end
            end

            %% Add the Significance bar
            sig_cell = {};
            sig_vec  = [];
            % for i = 1:num_vec
            %     for j = i+1:num_vec
            %         sig_cell{1,(i-1)*num_vec-i*(i-1)/2 + j - i} = [i,j];
            %         sig_vec = [sig_vec,ranksum(vec_mat(:,i),vec_mat(:,j))];
            %     end
            % end
            for i = 1:num_vec-1
                sig_cell{1,i} = [i,num_vec];
                sig_vec = [sig_vec,ranksum(vec_mat(:,i),vec_mat(:,num_vec))];
            end
            H = sigstar(sig_cell,sig_vec,0,scale_cmd);

            % keyboard

        case 'signed rank'
            boxplot(vec_mat,'Symbol','o','OutlierSize',3,'Colors',[0,0,0])
            ylim(Y_LIM)
            ax=gca;hold on;
            ax.LineWidth=1.1;
            ax.FontSize=22;
            ax.FontName='Arial';
            ax.FontWeight = 'bold';
            ax.XTickLabel=vec_name;
            ax.Title.String=title;
            ax.Title.FontSize=27;
            ax.YLabel.String=Y_label;
            if scale_cmd == 1
                ax.YScale = 'log';
            end

            %% Width of Line
            lineObj=findobj(gca,'Type','Line');
            for i=1:length(lineObj)
                lineObj(i).LineWidth=1;
                lineObj(i).MarkerFaceColor=[1,1,1].*.3;
                lineObj(i).MarkerEdgeColor=[1,1,1].*.3;
            end

            %% Coloring the box
            if ~isempty(colors)
                boxObj=findobj(gca,'Tag','Box');
                boxObj
                for i=1:length(boxObj)
                    patch(boxObj(i).XData,boxObj(i).YData,colors(i,:),'FaceAlpha',0.5,'LineWidth',1.1);
                end
            end

            %% Add the Significance bar
            sig_cell = {};
            sig_vec  = [];
            for i = 1:num_vec
                for j = i+1:num_vec
                    sig_cell{1,(i-1)*num_vec-i*(i-1)/2 + j - i} = [i,j];
                    sig_vec = [sig_vec,signrank(vec_mat(:,i)-vec_mat(:,j))];
                end
            end
            H = sigstar(sig_cell,sig_vec,0,scale_cmd);
            %% Add the paired line
            X=ones(size(vec_mat)).*(1:size(vec_mat,2));
            plot(X',vec_mat','Color',[0,0,0,.3],'Marker','o','MarkerFaceColor',[1,1,1].*.3,'MarkerEdgeColor',[1,1,1].*.3,'MarkerSize',3,'LineWidth',.6)
        case 'no_sig'
            boxplot(vec_mat,'Symbol','o','OutlierSize',3,'Colors',[0,0,0])
            ylim(Y_LIM)
            ax=gca;hold on;
            ax.LineWidth=1.1;
            ax.FontSize=22;
            ax.FontName='Arial';
            ax.FontWeight = 'bold';
            ax.XTickLabel=vec_name;
            ax.Title.String=title;
            ax.Title.FontSize=27;
            ax.YLabel.String=Y_label;
            if scale_cmd == 1
                ax.YScale = 'log';
            end

            %% Width of Line
            lineObj=findobj(gca,'Type','Line');
            for i=1:length(lineObj)
                lineObj(i).LineWidth=1;
                lineObj(i).MarkerFaceColor=[1,1,1].*.3;
                lineObj(i).MarkerEdgeColor=[1,1,1].*.3;
            end

            %% Coloring the box
            if ~isempty(colors)
                boxObj=findobj(gca,'Tag','Box');
                boxObj
                for i=1:length(boxObj)
                    patch(boxObj(i).XData,boxObj(i).YData,colors(i,:),'FaceAlpha',0.5,'LineWidth',1.1);
                end
            end


    end


end