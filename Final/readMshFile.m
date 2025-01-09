function [n_np, n_el, nodes_coords] = readMshFile(filename)
    % 打开文件
    fid = fopen(filename, 'r');
    if fid == -1
        error('无法打开文件: %s', filename);
    end
    
    % 初始化变量
    n_np = 0;
    n_el = 0;
    nodes_coords = [];
    
    % 逐行读取文件
    tline = fgetl(fid);
    while ischar(tline)
        % 查找 $Nodes
        if contains(tline, '$Nodes')
            % 读取下一行，获取节点数量 n_np
            tline = fgetl(fid);
            tokens = str2num(tline);  % 转化为数字
            if length(tokens) >= 2
                n_np = tokens(2);  % 第二个数字为节点数量
            end
            break;
        end
        tline = fgetl(fid);
    end
    
    % 读取节点坐标
    while ischar(tline)
        tline = fgetl(fid);
        tokens = str2num(tline);
        % 判断是否为节点坐标行（有 3 个数字）
        if length(tokens) == 3
            % 存储节点坐标
            nodes_coords = [nodes_coords; tokens(1), tokens(2)];  % 只存 x, y 坐标
        end
        
        % 如果遇到 $Elements，则停止读取节点
        if contains(tline, '$Elements')
            break;
        end
    end
    % 读取元素数量 n_el
    while ischar(tline)
        tline = fgetl(fid);
         if contains(tline, '$EndElements')
            break;
        end
        tokens = str2num(tline);
        % 判断是否为元素行（有 5 个数字）
        if length(tokens) == 5
            n_el = n_el + 1;  % 每有一行元素就加 1
        end
    end
    
    % 关闭文件
    fclose(fid);
end
