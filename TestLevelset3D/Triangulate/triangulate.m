
%[F1,V1] = geometry(LS);
[F,V,N] = triangulate(LS);

figure(5);
clf;
patch('Faces',F, 'Vertices',V, 'VertexNormals', N, 'FaceColor','red','EdgeColor','none');
view(3); axis vis3d;
camlight;
lighting gouraud;
