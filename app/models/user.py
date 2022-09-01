from flask_login import UserMixin

from app.db.db import get_db

class User(UserMixin):
    def __init__(self, id_, name, email):
        self.id = id_
        self.name = name
        self.email = email

    @staticmethod
    def get(user_id):
        """Retrieve user by id from database"""

        db = get_db()
        user = db.execute(
            "SELECT * FROM user WHERE id = ?", (user_id,)
        ).fetchone()
        if not user:
            return None

        user = User(
            id_=user[0], name=user[1], email=user[2]
        )
        return user

    @staticmethod
    def create(id_, name, email):
        """Add user to db"""

        db = get_db()
        if id and name and email:
            insert_fields = "id, name, email"
            insert_values = "?, ?, ?"
            insert_tuple = (id_, name, email)
        elif id and name:
            insert_fields = "id, name"
            insert_values = "?, ?"
            insert_tuple = (id_, name)
        elif id and email:
            insert_fields = "id, email"
            insert_values = "?, ?"
            insert_tuple = (id_, email)
        else:
            insert_fields = "id"
            insert_values = "?"
            insert_tuple = (id_,)
        db.execute(
            f"INSERT INTO user ({insert_fields}) "
            f"VALUES ({insert_values})",
            insert_tuple
        )
        db.commit()
