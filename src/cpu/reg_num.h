#define REG_NUM 32
#define ADDITIONAL_REGNUM 5

enum RenamedStatus {
    NOT_RENAMED = 0,
    BEING_RENAMED = 1,
    AFTER_RENAME = 2
}