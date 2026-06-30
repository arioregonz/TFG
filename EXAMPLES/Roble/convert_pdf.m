archivos = dir('*.fig');

for k = 1:length(archivos)

    nombre = archivos(k).name(1:end-4);  % quitar .fig

    openfig(archivos(k).name);

    print(gcf, nombre, '-dpng', '-r300');

    close;

end