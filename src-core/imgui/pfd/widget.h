#pragma once

#include "imgui/file_dialog/file_dialog.h"
#include "portable-file-dialogs.h"

class FileSelectWidget
{
public:
    FileSelectWidget(std::string label, std::string selection_text, bool directory = false, bool allow_url = false);
    ~FileSelectWidget();
    bool draw(std::string hint = "");
    bool isValid();
    std::string getPath();
    void setPath(std::string new_path);
    void setDefaultDir(std::string new_path);

private:
    std::string label, selection_text, id, btnid, default_dir, path;
    bool directory, waiting_for_res, allow_url;
    std::shared_ptr<satdump::FileDialog> fileselect;
    pfd::select_folder *dirselect;
    bool file_valid, url_valid;
};
