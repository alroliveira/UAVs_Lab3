function [path, nodes] = RRTStar_3D(q_start, q_goal, obs1, obs2, EPS, numNodes, bounds)
    q_start_node.coord = q_start;
    q_start_node.cost = 0;
    q_start_node.parent = 0;
    nodes(1) = q_start_node;

    for i = 1:numNodes
        q_rand = randInBounds(bounds);
        ndist = arrayfun(@(n) dist_3d(n.coord, q_rand), nodes);
        [val, idx] = min(ndist);
        q_near = nodes(idx);
        q_new.coord = steer3d(q_rand, q_near.coord, val, EPS);

        if noCollision(q_new.coord, q_near.coord, obs1) && noCollision(q_new.coord, q_near.coord, obs2)
            q_new.cost = dist_3d(q_new.coord, q_near.coord) + q_near.cost;

            q_nearest = [];
            r = 0.3;
            for j = 1:length(nodes)
                if noCollision(nodes(j).coord, q_new.coord, obs1) && noCollision(nodes(j).coord, q_new.coord, obs2) ...
                        && dist_3d(nodes(j).coord, q_new.coord) <= r
                    q_nearest = [q_nearest; nodes(j)];
                    line([q_near.coord(1), q_new.coord(1)], ...
                     [q_near.coord(2), q_new.coord(2)], ...
                     [q_near.coord(3), q_new.coord(3)], ...
                     'Color', 'k');  % Cinzento claro para os ramos "falhados"
                end
            end

            q_min = q_near;
            C_min = q_new.cost;
            for k = 1:length(q_nearest)
                cost_temp = q_nearest(k).cost + dist_3d(q_nearest(k).coord, q_new.coord);
                if cost_temp < C_min
                    q_min = q_nearest(k);
                    C_min = cost_temp;
                end
            end

            for j = 1:length(nodes)
                if all(nodes(j).coord == q_min.coord)
                    q_new.parent = j;
                end
            end

            nodes = [nodes q_new];

            if dist_3d(q_new.coord, q_goal) < EPS
                q_goal_node.coord = q_goal;
                q_goal_node.cost = q_new.cost + dist_3d(q_new.coord, q_goal);
                q_goal_node.parent = length(nodes);
                nodes = [nodes q_goal_node];
                break
            end
        end
    end

    % Reconstruir caminho
    path = [];
    q_end = nodes(end);
    while q_end.parent ~= 0
        path = [path; q_end.coord];
        q_end = nodes(q_end.parent);
    end
    path = [path; q_end.coord];
    path = flipud(path);
end
