function RenameSheet(FileName,SheetName)

%FileName：生成的excel表格的文件名称
%SheetName：生成的excel表格中的sheet表的名称

% 判断Excel是否已经打开，若已打开，就在打开的Excel中进行操作，否则就打开Excel
try
    % 若Excel服务器已经打开，返回其句柄Excel
    Excel = actxGetRunningServer('Excel.Application');
catch
    % 创建一个Microsoft Excel服务器，返回句柄Excel
    Excel = actxserver('Excel.Application');
end

% 设置Excel服务器为可见状态
Excel.Visible = 1;    % set(Excel, 'Visible', 1); 

% 若测试文件存在，打开该测试文件，否则，新建一个工作簿，并保存，文件名为测试.Excel
if exist(FileName,'file'); 
    Workbook = Excel.Workbooks.Open(FileName);
    % Workbook = invoke(Excel.Workbooks,'Open',filespec_user);
else
    Workbook = Excel.Workbooks.Add;
    % Workbook = invoke(Excel.Workbooks, 'Add');
    Workbook.SaveAs(FileName);
end

% 返回当前工作表句柄
Sheets = Excel.ActiveWorkbook.Sheets;    % Sheets = Workbook.Sheets;
Sheet1 = Sheets.Item(1);    % 返回第1个表格句柄
Sheet1.Activate;    % 激活第1个表格

Sheet1.Name = SheetName;

Workbook.Save   % 保存文档
Workbook.Close   % 保存文档

% 关闭Excel服务器
Excel.Visible = 0;    % set(Excel, 'Visible', 0); 
end