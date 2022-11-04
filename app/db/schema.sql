CREATE TABLE user (
  id TEXT PRIMARY KEY,
  name TEXT,
  email TEXT UNIQUE,
  is_guest INTEGER DEFAULT 0,
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- separator --
CREATE TABLE task (
  id TEXT NOT NULL,
  user_id TEXT NOT NULL,
  task_type TEXT NOT NULL,
  status TEXT NOT NULL,
  created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (id),
  FOREIGN KEY (user_id) REFERENCES user(id)
);

-- separator --
CREATE TRIGGER [UpdateLastTime]
    AFTER UPDATE
    ON task
    FOR EACH ROW
    WHEN NEW.updated_at < OLD.updated_at    --- this avoid infinite loop
BEGIN
    UPDATE task SET updated_at=CURRENT_TIMESTAMP WHERE id=OLD.id;
END;

-- separator -- 
INSERT INTO user (id, name, email, is_guest) VALUES('guest1', 'John', 'email@domain.com', 1);
