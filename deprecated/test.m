addpath('E:\Repository\Experiments\Lei.Zhang\utils\export_fig')
plot(cos(linspace(0, 7, 1000)));
set(gcf, 'Position', [100 100 150 150]);
saveas(gcf, 'test.png');
export_fig test2.png